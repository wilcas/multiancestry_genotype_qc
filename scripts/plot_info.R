library(data.table)
args <- commandArgs(trailingOnly=TRUE)
info <- lapply(
  dir(
    args[[1]],
    pattern=".*info.*.gz",
    full.names=TRUE
  ),
  fread
)

info <- rbindlist(info)
info_vec <- na.omit(as.numeric(info$Rsq))
breaks <- seq(0,1,0.001)
info_vec.cut <- cut(info_vec,breaks, right=F)
info_vec.freq <- table(info_vec.cut)
cum_freq <- c(cumsum(info_vec.freq))
png(args[[2]], type="cairo")
plot(breaks[-c(1)], cum_freq,type='l',xaxt='n', xlab = "info", ylab="cumulative frequency", main = 'Cumulative Frequency of R2 Scores')
axis(1, at = seq(0,1,0.01))
abline(v=0.05,lty='dashed')
abline(v=0.3,lty='dashed')
abline(v=0.95,lty='dashed')
dev.off()
