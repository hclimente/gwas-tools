#!/usr/bin/env Rscript
# RGWAS
# gwas.bed, gwas.bim, gwas.fam

library(snpStats)

load(${RGWAS})

write.plink(file.base = 'gwas',
            snps = gwas[['genotypes']],
            subject.data = gwas[['fam']],
            snp.data = gwas[['map']])