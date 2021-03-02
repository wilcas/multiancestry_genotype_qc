import pandas as pd
from itertools import combinations

rule pre_cleaning:
    input:
        bed="{dataset}.bed",
        bim="{dataset}.bim",
        fam="{dataset}.fam"

    output:
        bed="{dataset}.geno.bed",
        bim="{dataset}.geno.bim",
        fam="{dataset}.geno.fam"
    log:
        "{dataset}_precleaning.log"
    run:
        shell(
            """
            {config[plink]} --bfile {wildcards.dataset} --set-hh-missing --allow-no-sex --make-bed --geno 0.05 --out {wildcards.dataset}.geno > {log}
            """
        )


rule individual_qc:
    input:
        bed="{dataset}.geno.bed",
        bim="{dataset}.geno.bim",
        fam="{dataset}.geno.fam"
    output:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.fam"
    log:
        "{dataset}_individual_qc.log"
    run:
        shell(
            """
            cat {input.bim} | grep -E "A.*T|T.*A|G.*C|C.*G" > {wildcards.dataset}_tmp_ambiguous.txt
            cut -f 2 {input.bim} | sort | uniq -d > {wildcards.dataset}_dup_tmp.txt
            cat {wildcards.dataset}_tmp_ambiguous.txt {wildcards.dataset}_dup_tmp.txt | sort | uniq -d > {wildcards.dataset}_to_remove.txt
            echo "MAF" > {log}
            {config[plink]} --bfile "{wildcards.dataset}.geno" \
                --maf 0.01 \
                --exclude {wildcards.dataset}_to_remove.txt\
                --make-bed \
                --out "{wildcards.dataset}.geno.maf" | grep "removed" >> {log}

            #individual missingness
            echo "missingness" >> {log}
            {config[plink]} --bfile "{wildcards.dataset}.geno.maf" \
                --mind 0.02 \
                --make-bed \
                --out "{wildcards.dataset}.geno.maf.mind"  | grep "removed" >> {log}

            #checksex probes
            {config[plink]} --bfile "{wildcards.dataset}.geno.maf.mind" \
            --check-sex \
            --out "{wildcards.dataset}.geno.maf.mind" 

            # filter out samples that don't pass check
            cat {wildcards.dataset}.geno.maf.mind.sexcheck |
            grep -v "OK" | tail -n+2 > {wildcards.dataset}.sex_failed.txt
            
            echo "sex" >> {log}
            {config[plink]} --bfile "{wildcards.dataset}.geno.maf.mind" \
                --remove {wildcards.dataset}.sex_failed.txt \
                --make-bed \
                --out "{wildcards.dataset}.geno.maf.mind.sex_check" 

            # Compute method of moments heritability F
            {config[plink]} --bfile "{wildcards.dataset}.geno.maf.mind.sex_check" \
                --het \
                --out "{wildcards.dataset}.geno.maf.mind.sex_check"
            awk   "{{ if (\$6 > 0.2) {{print \$1,\$2 }} }}" {wildcards.dataset}.geno.maf.mind.sex_check.het > {wildcards.dataset}.f_failed.txt
            echo "Fhet" >> {log}
            {config[plink]} --bfile "{wildcards.dataset}.geno.maf.mind.sex_check"\
                --remove {wildcards.dataset}.f_failed.txt \
                --make-bed \
                --out "{wildcards.dataset}.geno.maf.mind.sex_check.het_filter" 
            """
        )


rule ibd_calculation:
    input:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.fam",
        clst="{dataset}.clst"
    output:
        gen="{dataset}.genome.gz",
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.fam"
    run:
        shell(
            """
            {config[plink]} --bfile {wildcards.dataset}.geno.maf.mind.sex_check.het_filter\
                --within "{input.clst}" \
                --maf 0.05 \
                --indep-pairwise 200 100 0.2 \
                --make-bed \
                --out {wildcards.dataset}.geno.maf.mind.sex_check.het_filter.ld1.maf05
            """
        )
        clst = pd.read_csv(input.clst, sep="\s+", header=None)
        pairs = [map(str, comb) for comb in combinations(clst.loc[:,2].unique().tolist(),2)]
        for (p1,p2) in pairs:
            shell(
                """
                {config[plink]}  --bfile {wildcards.dataset}.geno.maf.mind.sex_check.het_filter.ld1.maf05\
                    --within {input.clst}\
                    --keep-cluster-names "{p1}" "{p2}"\
                    --fst\
                    --out {wildcards.dataset}_{p1}_{p2}
                awk "{{if (\$5 > 0.2 || \$5 < -0.2) {{print \$2}}}}" {wildcards.dataset}_{p1}_{p2}.fst > {wildcards.dataset}_{p1}_{p2}_remove.txt
                """
            )
        shell("cat {wildcards.dataset}*_remove.txt | sort | uniq > {wildcards.dataset}_fst_remove.txt")
        shell(
            """
            {config[plink]} --bfile {wildcards.dataset}.geno.maf.mind.sex_check.het_filter.ld1.maf05\
                --exclude {wildcards.dataset}_fst_remove.txt\
                --genome gz full unbounded nudge\
                --out {wildcards.dataset}
            """)
        genome = pd.read_csv(wildcards.dataset + ".genome.gz", sep="\s+")
        genome.loc[genome.PI_HAT > 0.2,["FID1","IID1"]].to_csv(wildcards.dataset + "_ibd_remove.txt",index=False,sep="\t")
        shell(
            """
            {config[plink]} --bfile {wildcards.dataset}.geno.maf.mind.sex_check.het_filter\
                --remove {wildcards.dataset}_ibd_remove.txt\
                --make-bed\
                --out {wildcards.dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter
            """
        )

rule pre_ancestry_pca:
    input:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.fam"
    output:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.fam",
        evec="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.eigenvec",
        evalu="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.eigenval"
    run:
        start=wildcards.dataset + ".geno.maf.mind.sex_check.het_filter.ibd_filter"
        cur_geno=wildcards.dataset + ".geno.maf.mind.sex_check.het_filter.ibd_filter"
        # Filter unreliable regions
        shell(
            """
            {config[plink]} --bfile "{cur_geno}" \
            --maf 0.05 \
            --autosome \
            --make-bed \
            --out "{cur_geno}.autosome"
            """
        )
        cur_geno+=".autosome"
        shell("cat {config[LD]} {config[mhc]} > {wildcards.dataset}.pca_to_remove.txt")
        shell(
            """
            {config[plink]} --bfile "{cur_geno}" \
                --exclude 'range' {wildcards.dataset}.pca_to_remove.txt \
                --make-bed \
                --out "{cur_geno}.exclude_regions"
            """
        )
        cur_geno+=".exclude_regions"
        # remove remaining regions in LD
        shell(
            """
            {config[plink]} --bfile "{cur_geno}" \
                --indep-pairwise 200 100 0.2 \
                --make-bed \
                --out "{cur_geno}.LD_1"

            {config[plink]} --bfile "{cur_geno}.LD_1" \
                --indep-pairwise 200 100 0.2 \
                --make-bed \
                --out "{start}.pre_pca"
            {config[plink]} --bfile "{start}.pre_pca" \
                --pca\
                --out "{start}"
            """
        )


rule pc_project:
    input:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.fam"
    output:
        kg_pcs="{dataset}.1kg_pca.pcs",
        proj="{dataset}.kg_projection.txt"
        
    run:
        cur_geno=wildcards.dataset + ".geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca"
        # match data to 1000 Genomes based on position
        shell("{config[Rscript]} scripts/merge_bim.R {config[kg_data_pre_pca]}.bim {input.bim} {wildcards.dataset}_filter.txt {wildcards.dataset}_matched.bim")
        matched_data = wildcards.dataset + "_matched.bim"
        filtered = wildcards.dataset + "_filter.txt"
        shell(
            """
            {config[plink]} --bfile {cur_geno} \
                --extract {filtered} \
                --make-bed \
                --out {wildcards.dataset}.kg_matched

            {config[plink]} --bed {wildcards.dataset}.kg_matched.bed --bim {matched_data} --fam {wildcards.dataset}.kg_matched.fam \
                --make-bed \
                --out {wildcards.dataset}.kg_matched.unified
            """
        )
          
        cur_geno = wildcards.dataset + ".kg_matched.unified"
        # filter 1000 genomes data to snps remaining in dataset
        shell(
            """
            {config[plink]} --bfile {config[kg_data_pre_pca]}\
                --extract {cur_geno}.bim\
                --make-bed\
                --out {wildcards.dataset}.1kg_pca
            """
        )
        kg_data=wildcards.dataset + ".1kg_pca"
        # PCA
        shell(
            """
            {config[flashpca]} --bfile {kg_data} \
                --outmeansd {kg_data}.meansd\
                --outload {kg_data}.loadings\
                --outpc {output.kg_pcs}
            {config[flashpca]} --project\
                --bfile {cur_geno}\
                --inmeansd {kg_data}.meansd\
                --inload {kg_data}.loadings\
                --outproj {output.proj}
            """
        )

rule ancestry_assignment:
    input:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.fam",
        clst="{dataset}.clst",
        kg_pcs="{dataset}.1kg_pca.pcs",
        proj="{dataset}.kg_projection.txt",
        ancestry_clusters = "{dataset}.ancestry_clusters.txt"
    output:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.fam",
        clst="{dataset}_assigned_ancestry.clst"
    run:
        # compute mahalanobis distances
        shell("{config[Rscript]} scripts/mahalanobis.R {input.kg_pcs} {input.proj} {input.ancestry_clusters} {config[kg_panel]} {output.clst}")
        cur_geno = input.bed[0:-4]
        shell(
            """
            {config[plink]} --bfile {cur_geno} \
                --keep {output.clst}\
                --make-bed\
                --out {cur_geno}.ancestry
            """
        )


rule batch_checking:
    input:
        batches="{dataset}_sample_batches.txt",
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.fam",
        clst="{dataset}_assigned_ancestry.clst",
        pca_bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bed",
        pca_bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bim",
        pca_fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.fam"
    output:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.fam",
        pre_impute_pcs="{dataset}.preimpute.pcs"
    run:
        # construct covariate file
        cur_geno = input.bed[0:-4]
        pca_bfile = input.pca_bed[0:-4] + ".ancestry"
        shell(
            """
            {config[plink]} --bfile {cur_geno}\
                --keep {input.clst}\
                --make-bed\
                --out {pca_bfile}

            {config[flashpca]} --bfile {pca_bfile} \
                --outpc {wildcards.dataset}.preimpute.pcs
            """
        )
        cur_fam = pd.read_csv(input.fam,sep="\s+",header=None)
        pcs = pd.read_csv(wildcards.dataset + ".preimpute.pcs", sep="\s+")
        covar = pd.concat([pcs,cur_fam[[4]]], axis=1)
        tmp_col = covar.columns.values
        tmp_col[-1] = "sex"
        covar.columns = tmp_col
        covar_file = wildcards.dataset+ "_batch_covar.cov"
        covar.to_csv(covar_file, index=False, sep= " ")
        # construct batch fam files
        batch_manifest = pd.read_csv(input.batches,sep="\s+")
        if len(batch_manifest.index) <= 1:
            shell(
                """
                cp {input.bed} {output.bed}
                cp {input.bim} {output.bim}
                cp {input.fam} {output.fam}
                """
            )
        else:
            for col in batch_manifest:
                for val in batch_manifest[col].unique():
                    cur_fam[cur_fam.columns[5]] = [1 if x == val else 2 for x in batch_manifest[col]]
                    cur_fam_file = f"{wildcards.dataset}_{str(val)}_batch.fam" 
                    cur_fam.to_csv(cur_fam_file, index=False, sep = " ",header=False)
                    shell(
                        """
                        {config[plink]} --bed {input.bed}\
                            --bim {input.bim}\
                            --fam {cur_fam_file}\
                            --within {input.clst}\
                            --logistic \
                            --covar {covar_file} \
                            --out {wildcards.dataset}_{val}
                        """
                    )
                    # Select snps significantly associated with batch
                    result = pd.read_csv(wildcards.dataset + "_" + str(val) + ".assoc.logistic", sep="\s+")
                    # remove snps significantly associated with batch
                    sig_assoc = result.loc[(result.TEST == "ADD") , :]
                    sig_assoc = sig_assoc.loc[sig_assoc.P < (0.05 / len(sig_assoc.index)), "SNP"]
                    sig_assoc.to_csv(wildcards.dataset + "_" + str(val) + "_sig_assoc.txt")
            # filter results
            shell(
                """
                cat {wildcards.dataset}*sig_assoc.txt > {wildcards.dataset}_all_batch_sig.txt
                {config[plink]} --bfile {cur_geno} \
                    --exclude {wildcards.dataset}_all_batch_sig.txt\
                    --make-bed \
                    --out {cur_geno}.batch
                """
            )
        
        
rule imputation_batch_dosages:
    input:
        gen="{data_dir}/raw_data.imputed.r2_30.maf_mismatch.gen",
        sample="{data_dir}/raw_data.imputed.r2_30.maf_mismatch.sample",
        batches="{data_dir}/imputed_sample_batches.txt",
        pca_bed="{data_dir}/raw_data.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bed",
        pca_bim="{data_dir}/raw_data.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.bim",
        pca_fam="{data_dir}/raw_data.geno.maf.mind.sex_check.het_filter.ibd_filter.pre_pca.fam"
    output:
        gen="{data_dir}/raw_data.imputed.r2_30.maf_mismatch.batch.gen",
        sample="{data_dir}/raw_data.imputed.r2_30.maf_mismatch.batch.sample"
    run:
        pca_bfile = input.pca_bed[0:-4]
        # construct genotyping PCs
        shell(
            """
            {config[flashpca]} --bfile {pca_bfile} \
                --outpc {wildcards.data_dir}/raw_data.preimpute.pcs
            """
        )
        # load covariates
        cur_fam = pd.read_csv(input.pca_fam,sep="\s+",header=None)
        pcs = pd.read_csv(wildcards.data_dir + "/raw_data.preimpute.pcs", sep="\s+")
        covar = pd.concat([pcs,cur_fam[[4]]], axis=1)
        tmp_col = covar.columns.values
        tmp_col[-1] = "sex"
        covar.columns = tmp_col
        covar_file = wildcards.data_dir+"/raw_data_batch_covar.cov"
        covar.to_csv(covar_file, index=False, sep= " ")
        # construct batch sample files
        batch_manifest = pd.read_csv(input.batches,sep="\s+")
        cur_sample = pd.read_csv(input.sample,sep="\s+")
        for col in batch_manifest:
            for val in batch_manifest[col].unique():
                sex = ["D"]
                sex.extend(cur_fam[4].tolist())
                cur_sample["sex"] = sex
                tmp = ["B"]
                tmp.extend([1 if x == val else 0 for x in batch_manifest[col]])
                cur_sample[col + str(val)] = tmp
                cur_sample_file = f"{wildcards.data_dir}/raw_data_{str(val)}_batch.sample" 
                cur_sample.to_csv(cur_sample_file, index=False, sep = " ")
                shell(
                    """
                    {config[plink2]} --gen {input.gen} ref-unknown\
                        --sample {cur_sample_file}\
                        --glm\
                        --covar {covar_file} \
                        --out {wildcards.data_dir}/imputed_dosage_batch_{val}
                    """
                )
                # Select snps significantly associated with batch
                sig_assoc = pd.read_csv(wildcards.data_dir + "/imputed_dosage_batch_" + str(val) + ".glm.logistic.hybrid", sep="\s+")
                # remove snps significantly associated with batch
                sig_assoc = sig_assoc.loc[sig_assoc.P < (0.05 / len(sig_assoc.index)), "ID"]
                sig_assoc.to_csv(wildcards.data_dir + "/raw_data_" + str(val) + "_sig_assoc.txt")
        # filter results
        shell(
            """
            cat {wildcards.data_dir}/*sig_assoc.txt > {wildcards.data_dir}/all_batch_sig.txt
            {config[plink2]} --gen {input.gen} ref-unknown\
                --sample {input.sample}\
                --exclude {wildcards.data_dir}/all_batch_sig.txt\
                --recode oxford\
                --out {data_dir}/raw_data.imputed.r2_30.maf_mismatch.batch
            """
        )
#rule imputation_batch_hard_calls:
#    input:
#        "{dataset}_"
#        batches="{dataset}_imputed_sample_batches.txt"
#    output:
#    
#    run:
    
rule hardy_weinberg:
    input:
        clst="{dataset}_assigned_ancestry.clst",
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.fam"
    output:
        "{dataset}_all_hwe.hwe.gz",
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.hwe.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.hwe.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.hwe.fam"
    run:
        cur_geno = input.bed[0:-4]
        ancestry = pd.read_csv(input.clst, sep="\s+", header=None)
        for clst in ancestry[2].unique():
            filter_file = wildcards.dataset +"_"+str(clst)+".txt"
            ancestry[ancestry[2] == clst].to_csv(filter_file,index=False,header=False,sep=" ")
            shell(
                """
                {config[plink]} --bfile {cur_geno}\
                    --keep {filter_file}\
                    --hardy gz\
                    --out {filter_file}
                """
            )
        #remove hardy_weinberg p
        shell("cat {wildcards.dataset}*.hwe.gz > {wildcards.dataset}_all_hwe.hwe.gz")
        shell("gunzip {wildcards.dataset}_all_hwe.hwe.gz")
        shell("head -n 1 {wildcards.dataset}_all_hwe.hwe > {wildcards.dataset}.tmp")
        shell("cat {wildcards.dataset}_all_hwe.hwe | grep -v 'O(HET)' >> {wildcards.dataset}.tmp")
        shell("mv {wildcards.dataset}.tmp {wildcards.dataset}_all_hwe.hwe; gzip {wildcards.dataset}_all_hwe.hwe")
        hw_all = pd.read_csv(wildcards.dataset + "_all_hwe.hwe.gz", sep="\s+")
        hw_all.fillna(0, inplace=True)
        print(hw_all.head)
        hw_all[hw_all.P < 1e-10]["SNP"].to_csv(wildcards.dataset + "_rm_hwe.txt", header=False, index=False, sep=" ")
        shell(
            """
            {config[plink]} --bfile {cur_geno}\
                --exclude {wildcards.dataset}_rm_hwe.txt\
                --make-bed\
                --out {cur_geno}.hwe
            """
        )

rule align_to_hg19:
    input:
        bed="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.hwe.bed",
        bim="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.hwe.bim",
        fam="{dataset}.geno.maf.mind.sex_check.het_filter.ibd_filter.ancestry.batch.hwe.fam"
    output:
        "{dataset}.pre_imputation.vcf.gz",
        "{dataset}.pre_imputation.vcf.gz.csi",
        expand("{{dataset}}.pre_imputation_chr{i}.vcf.gz", i=[str(x) for x in range(1,23)]+["X"],allow_missing=True)
    run:
        cur_geno = input.bed[0:-4]
        shell(
            """
            {config[plink]} --bfile {cur_geno}\
                --snps-only just-acgt\
                --set-hh-missing\
                --output-chr MT\
                --recode vcf-iid\
                --out {cur_geno}
            {config[bcftools]} +fixref {cur_geno}.vcf -Oz -o {wildcards.dataset}.pre_imputation.vcf.gz -- -f {config[hg19]} -m top
            {config[bcftools]} index {wildcards.dataset}.pre_imputation.vcf.gz
            """
        )
        for i in [str(x) for x in range(1,23)]+["X"]:
            shell("{config[bcftools]} view {wildcards.dataset}.pre_imputation.vcf.gz -r {i} -Oz -o {wildcards.dataset}.pre_imputation_chr{i}.vcf.gz")

# threshold for MAF difference between our study and a reference population
# and MAF  <0.01 threshold
rule merge_imputed_and_filter:
    input:
        dosage=expand("{{data_dir}}/chr{i}.dose.vcf.gz", i= [str(x) for x in range(1,23)]+["X"],allow_missing=True),
        info=expand("{{data_dir}}/chr{i}.info.gz", i= [str(x) for x in range(1,23)]+["X"],allow_missing=True)
    output:
        #"{data_dir}/all_imputed_unfilter_r2_{thresh}.vcf.gz",
        #"{data_dir}/all_imputed_r2_{thresh}.vcf.gz",
        "{data_dir}/all_imputed_r2_{thresh}_rsid.vcf.gz"#,
        #"{data_dir}/info_plot_r2_{thresh}.png"
    run:
        #plot info scores
        data_dir=wildcards.data_dir
        thresh= "0." + wildcards.thresh
        #shell("{config[Rscript]} scripts/plot_info.R {data_dir} {data_dir}/info_plot_r2_{wildcards.thresh}.png")
        shell(
            """
            cd {data_dir}
            for full_f in {input.dosage}; do
                f=$(basename $full_f)
                (
                    {config[bcftools]} filter -e "R2<{thresh}" --threads {config[threads]} -O z -o ${{f%%.*}}_r2_{wildcards.thresh}.vcf.gz $f 
                    {config[tabix]} -p vcf ${{f%%.*}}_r2_{wildcards.thresh}.vcf.gz 
                    {config[bcftools]} annotate -a {config[kg_vars]} -c CHROM,ID -O z --threads {config[threads]} -o ${{f%%.*}}_r2_{wildcards.thresh}_rsid.vcf.gz ${{f%%.*}}_r2_{wildcards.thresh}.vcf.gz 
                ) &
            done 
            wait
            cd ..
            {config[bcftools]} concat  {data_dir}/chr{{1..22}}_r2_{wildcards.thresh}_rsid.vcf.gz {data_dir}/chrX_r2_{wildcards.thresh}_rsid.vcf.gz -o {data_dir}/all_imputed_r2_{wildcards.thresh}_rsid.vcf.gz
            """
        )

rule vcf_to_oxford:
    input:
        vcf="{data_dir}/all_imputed_r2_{thresh}_rsid.vcf.gz"
    output:
        gen="{data_dir}/all_imputed_r2_{thresh}_rsid.gen",
        sample="{data_dir}/all_imputed_r2_{thresh}_rsid.sample"
    run:
        shell(
            """
            {config[plink2]} --vcf {input.vcf} dosage=HDS\
                --double-id \
                --max-alleles 2\
                --snps-only just-acgt\
                --recode oxford\
                --out {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid
            """
        )

rule dosage_maf:
    input:
        gen="{data_dir}/all_imputed_r2_{thresh}_rsid.gen",
        sample="{data_dir}/all_imputed_r2_{thresh}_rsid.sample",
        clst="{data_dir}/raw_data_assigned_ancestry.clst",
        ancestry_clusters="{data_dir}/raw_data.ancestry_clusters.txt"
    output:
        "{data_dir}/raw_data.imputed.r2_{thresh}.maf_mismatch.traw",
        "{data_dir}/raw_data.imputed.r2_{thresh}.maf_mismatch.gen",
        "{data_dir}/raw_data.imputed.r2_{thresh}.maf_mismatch.sample"
    run:
        ancestry=pd.read_csv(input.clst,sep="\s+",header=None)
        ancestry[0] = ancestry[1]
        for pop in ancestry[2].unique(): 
            panel = pd.read_csv(config['kg_panel'],sep="\s+")
            kg_filter = wildcards.data_dir + "/kg_"+pop+"_impute_filter.txt"
            panel = panel[(panel['pop'] == pop) | (panel['super_pop'] == pop)]
            if panel.shape[0] > 1:
                panel.insert(0,0,0)
                panel.to_csv(kg_filter, sep=" ", header=False, index=False)
                shell(
                    """
                    {config[plink2]} --bfile {config[kg_data]}\
                        --freq \
                        --keep {kg_filter}\
                        --out {wildcards.data_dir}/imputed_freq_dosages_1kg_{pop}
                    """
                )
                # Compute MAF in my data
                pop_filter = wildcards.data_dir +"/"+pop+"_impute_filter.txt"
                ancestry[ancestry[2] == pop].to_csv(pop_filter,sep=" ",header=False, index=False)
                shell(
                    """
                    {config[plink2]} --data {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid ref-unknown\
                        --freq \
                        --keep {pop_filter}\
                        --out {wildcards.data_dir}/imputed_freq_dosages_{pop}
                    """
                )
                # Compare frequencies and store
                kg_freq = pd.read_csv(wildcards.data_dir +"/imputed_freq_dosages_1kg_"+pop+".afreq",sep="\s+")
                my_freq = pd.read_csv(wildcards.data_dir+"/imputed_freq_dosages_"+pop+".afreq",sep="\s+")
                merged = my_freq.merge(kg_freq, on=['#CHROM','ID'])
                merged = merged[merged["ALT_FREQS_x"] != 0]
                merged = merged[(abs(merged["ALT_FREQS_x"] - merged["ALT_FREQS_y"]) > 0.15) | (merged["ALT_FREQS_x"] < 0.01)]
                merged.to_csv(wildcards.data_dir+"/imputed_freq_merged_"+pop+".txt",sep="\t",index=False)
                to_remove = wildcards.data_dir+"/to_filter_freq"+pop+".txt"
                merged["ID"].to_csv(to_remove, sep=" ", index=False,header=False)
            outfile =output[0][:-5]
        shell(
            """
            cat {wildcards.data_dir}/to_filter_freq*.txt > {wildcards.data_dir}/all_freq_to_remove.txt
            {config[plink2]} --data {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid ref-unknown\
                --exclude {wildcards.data_dir}/all_freq_to_remove.txt\
                --recode A-transpose\
                --out {outfile} 
              
            {config[plink2]} --data {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid ref-unknown\
                --exclude {wildcards.data_dir}/all_freq_to_remove.txt\
                --recode oxford\
                --out {outfile}
            """
        )


            
        
# Calculate frequencies of reference and filter variants with MAF different in my data
rule hard_call_maf:
    input:
        gen="{data_dir}/all_imputed_r2_{thresh}_rsid.gen",
        sample="{data_dir}/all_imputed_r2_{thresh}_rsid.sample",
        clst="{data_dir}/raw_data_assigned_ancestry.clst",
        ancestry_clusters="{data_dir}/raw_data.ancestry_clusters.txt"
    output:
        "{data_dir}/raw_data.imputed.r2_{thresh}.hard_call.maf_mismatch.bed",
        "{data_dir}/raw_data.imputed.r2_{thresh}.hard_call.maf_mismatch.bim",
        "{data_dir}/raw_data.imputed.r2_{thresh}.hard_call.maf_mismatch.fam"
        
    run:
        # hard call
        shell(
            """
            {config[plink2]} --data {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid ref-unknown\
                --hard-call-threshold 0.1\
                --make-bed\
                --out {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid_hard_call
            """
        )
        ancestry=pd.read_csv(input.clst,sep="\s+",header=None)
        for pop in ancestry[2].unique(): 
            panel = pd.read_csv(config['kg_panel'],sep="\s+")
            kg_filter = wildcards.data_dir + "/kg_"+pop+"_impute_filter.txt"
            panel = panel[(panel['pop'] == pop) | (panel['super_pop'] == pop)]
            panel.insert(0,0,0)
            panel.to_csv(kg_filter, sep=" ", header=False, index=False)
            shell(
                """
                {config[plink2]} --bfile {config[kg_data]}\
                    --freq \
                    --keep {kg_filter}\
                    --out {wildcards.data_dir}/imputed_freq_dosages_1kg_{pop}
                """
            )
            # Compute MAF in my data
            pop_filter = wildcards.data_dir +"/"+pop+"_impute_filter.txt"
            ancestry[ancestry[2] == pop].to_csv(pop_filter,sep=" ",header=False, index=False)
            shell(
                """
                {config[plink2]} --bfile {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid_hard_call\
                    --freq \
                    --keep {pop_filter}\
                    --out {wildcards.data_dir}/imputed_freq_dosages_hard_call_{pop}
                """
            )
            # Compare frequencies and store
            kg_freq = pd.read_csv(wildcards.data_dir +"/imputed_freq_dosages_1kg_"+pop+".afreq",sep="\s+")
            my_freq = pd.read_csv(wildcards.data_dir+"/imputed_freq_dosages_hard_call_"+pop+".afreq",sep="\s+")
            merged = my_freq.merge(kg_freq, on=['#CHROM','ID'])
            merged = merged[merged["ALT_FREQS_x"] != 0]
            merged = merged[(abs(merged["ALT_FREQS_x"] - merged["ALT_FREQS_y"]) > 0.15) | (merged["ALT_FREQS_x"] < 0.01)]
            merged.to_csv(wildcards.data_dir+"/imputed_freq_merged_"+pop+".txt",sep="\t",index=False)
            to_remove = wildcards.data_dir+"/to_filter_freq"+pop+".txt"
            merged["ID"].to_csv(to_remove, sep=" ", index=False,header=False)
            outfile =output[0][0:-4]
        shell(
            """
            cat {wildcards.data_dir}/to_filter_freq*.txt > {wildcards.data_dir}/all_freq_to_remove.txt
            {config[plink2]} --bfile {wildcards.data_dir}/all_imputed_r2_{wildcards.thresh}_rsid_hard_call\
                --exclude {wildcards.data_dir}/all_freq_to_remove.txt\
                --make-bed\
                --out {outfile}
            """
        )
