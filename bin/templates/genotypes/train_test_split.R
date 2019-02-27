#!/usr/bin/env Rscript
# RGWAS, K, I
# split_*.RData

load('${RGWAS}')
set.seed(0)

# permute samples
n <- nrow(gwas[['genotypes']])
perm <- sample(n)
  
gwas[['genotypes']] <- gwas[['genotypes']][perm,]
gwas[['fam']] <- gwas[['fam']][perm,]

# split folds  
K <- cut(seq(1, n), breaks = ${K}, labels = FALSE)

# train
train <- list()
train[['genotypes']] <- gwas[['genotypes']][K != ${I}, ]
train[['fam']] <- gwas[['fam']][K != ${I}, ]
train[['map']] <- gwas[['map']]
save(train, file = paste0('train_', ${I}, '.RData'))

# test
test <- list()
test[['genotypes']] <- gwas[['genotypes']][K != ${I}, ]
test[['fam']] <- gwas[['fam']][K != ${I}, ]
test[['map']] <- gwas[['map']]
save(test, file = paste0('test_', ${I}, '.RData'))