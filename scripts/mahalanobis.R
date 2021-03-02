library(tidyverse)
# Load resources
argv <- commandArgs(trailingOnly=TRUE)
kg_pcs <- read.delim(argv[[1]],sep="",stringsAsFactors = F)
data_projected <- read.delim(argv[[2]],sep="",stringsAsFactors = F)
clusters <- read.delim(argv[[3]],sep="",stringsAsFactors = F) # population is_super_pop
kg_panel <- read.delim(argv[[4]], sep="",stringsAsFactors = F)
out_name <- argv[[5]]

#get sample indices for selected clusters
samples <- lapply(
  1:nrow(clusters),
  function(i){
    if(clusters[i,2]){
      kg_panel$sample[kg_panel$super_pop == clusters[i,1]]
    }else{
      kg_panel$sample[kg_panel$pop == clusters[i,1]]
    }
  }
)

#filter outlier samples
pop_samples <- list()
for(i in 1:length(clusters$population)){
  pop_samples[[clusters$population[i]]] <-  kg_pcs %>%
    filter(IID %in% samples[[i]]) %>% 
    select(contains("PC")) %>%
    filter(sqrt(mahalanobis(.,colMeans(.),cov(.))) < 4)
}

#classify my samples and identify "ancestral" outliers
just_pcs <- data_projected %>% select(contains("PC"))
pop_dists <- lapply(
  names(pop_samples),
  function(pop){
    lab <- paste0(pop,"_dist")
    res <- data.frame(tmp=
      sqrt(
        mahalanobis(
          just_pcs,
          colMeans(pop_samples[[pop]]),
          cov(pop_samples[[pop]])
        )
      )
    )
    colnames(res)[1] <- lab
    return(res)
  }
)
pop_dists <- cbind(data_projected[,c("FID","IID")],do.call(cbind,pop_dists))
print(pop_dists)
assign_cluster <- pop_dists %>%
  mutate_at(vars(contains("_dist")),list("sd"=~ . > (4 * .))) %>%
  filter_at(vars(contains("sd")),all_vars(!.))%>%
  select(-contains("sd")) %>% 
  pivot_longer(contains("dist")) %>%
  mutate(name=gsub("_dist","",name)) %>%
  group_by(FID,IID) %>%
  summarize(assigned=name[which.min(value)])
reclassified_data <- data_projected %>% 
  left_join(assign_cluster,by=c("FID","IID")) %>%
  replace_na(list(assigned="excluded"))
write.table(
  reclassified_data %>%
    select(FID,IID,assigned) %>%
    filter(assigned != "excluded"),
  out_name,
  sep=" ",
  col.names = F,
  row.names = F,
  quote = F
)
                   

