#!/usr/bin/env Rscript
#BED, BIM, FAM
# x.npy, y.npy

library(snpStats)
library(RcppCNPy)

gwas <- read.plink("${BED}", "${BIM}", "${FAM}")

X <- as(gwas\$genotypes, "numeric")
Y <- gwas\$fam\$affected

npySave('x.npy', X)
npySave('y.npy', Y)