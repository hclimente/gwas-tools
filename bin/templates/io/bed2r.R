#!/usr/bin/env Rscript
# BED, BIM, FAM
# gwas.RData

library(snpStats)

gwas <- read.plink("${BED}", "${BIM}", "${FAM}")

save(gwas, file = 'gwas.RData')