#!/usr/bin/env Rscript
# PED, MAP
# out.sample, out.gen
  
library(snpStats)

gwas <- read.pedfile("${PED}", snps = "${MAP}")

save(gwas, file = 'gwas.RData')