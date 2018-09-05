#!/usr/bin/Rscript
# Read newline-separated numbers, calculate mean + 95% confidence intervals using basic bootstrap

library(boot)
set.seed(1)

#a <- c(rep(c(0,1), times=c(1000,1000)))
a <- read.table(file('stdin'), header=FALSE, sep="\t", row.names=NULL)$V1

bootmean <- boot(data=a, statistic=function(x, i) {mean(x[i])}, R=1000)

bootci <- boot.ci(boot.out = bootmean, type = c("basic"))
lo <- bootci$basic[4][1]
hi <- bootci$basic[5][1]
#bootci <- boot.ci(boot.out = bootmean, type = c("perc"))
#lo <- bootci$perc[4][1]
#hi <- bootci$perc[5][1]

cat(bootmean$t0, lo, hi, sep='\n')
