#!/usr/bin/env Rscript
library(martini)
library(igraph)
library(tidyverse)

load("$RGWAS")
load("$RNET")

params <- capture.output(
  subnet <- scones.cv(gwas, net,
                      score = "${SCORE}",
                      criterion = "${CRITERION}"),
  type = "message"
) %>% tail(n = 2) %>% lapply(strsplit, '=') %>% unlist %>% .[c(F,T)] %>% as.numeric() %>% log10

etas <- 10^seq(params[1] - 1.5, params[1] + 1.5, length.out = 10)
lambdas <- 10^seq(params[2] - 1.5, params[2] + 1.5, length.out = 10)
subnet <- scones.cv(gwas, 
                    net,
                    score = "${SCORE}",
                    criterion = "${CRITERION}",
                    etas = etas,
                    lambdas = lambdas)

cones <- martini:::get_snp_modules(gwas, subnet)
cones[['c']] <- martini:::snp_test(gwas, covars, "${SCORE}")
cones[['selected']] <- cones[['snp']] %in% names(V(subnet))

write_tsv(cones, 'cones.tsv')
