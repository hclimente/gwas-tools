#!/usr/bin/env Rscript
library(martini)
library(tidyverse)

load("${RGWAS}")
load("${RNET}")

# make a exploratory run to get the best parameters
params <- capture.output(
    cones <- search_cones(gwas, net, associationScore = "${SNP_SCORE}",
                          modelScore = "${MODEL_SCORE}", encoding = "${ENCODING}")
                        ) %>% 
        tail(n = 2) %>% 
        lapply(strsplit, ' = ') %>% 
        unlist %>% .[c(F,T)] %>% 
        as.numeric() %>% 
        log10

# optimize the parameters
etas <- 10^seq(params[1] - 1, params[1] + 1, length.out = 10)
lambdas <- 10^seq(params[2] - 1, params[2] + 1, length.out = 10)

cones <- search_cones(gwas, net,
                      associationScore = "${SNP_SCORE}",
                      modelScore = "${MODEL_SCORE}",
                      encoding = "${ENCODING}",
                      etas = etas,
                      lambdas = lambdas)
  
write_tsv(cones, 'cones.tsv')