require(data.table)
args <- commandArgs(trailingOnly = TRUE)
reference <- fread(args[[1]])
to_align <- fread(args[[2]])
filter_file <- args[[3]]
out_file <- args[[4]]
merged <- merge(to_align,reference,by=c("V1","V3","V4"))
fwrite(merged[,.(V2.x)],filter_file,quote = F,row.names=F,col.names = F,sep = " ")
fwrite(merged[,.(V1,V2.y,V3,V4,V5.y,V6.y)],out_file,quote = F,row.names=F,col.names = F,sep = " ")