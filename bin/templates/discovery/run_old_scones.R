#!/usr/bin/env Rscript
library(martini)
library(igraph)
library(tidyverse)

load("$RGWAS")
load("$RNET")

params <- capture.output(
  cones <- search_cones(gwas, net,
                        associationScore = "${SCORE}",
                        modelScore = "${}CRITERION}")
) %>% tail(n = 2) %>% lapply(strsplit, ' = ') %>% unlist %>% .[c(F,T)] %>% as.numeric() %>% log10

etas <- 10^seq(params[1] - 1.5, params[1] + 1.5, length.out = 10)
lambdas <- 10^seq(params[2] - 1.5, params[2] + 1.5, length.out = 10)
cones <- search_cones(gwas, net,
                      associationScore = "${SCORE}",
                      modelScore = "${CRITERION}",
                      etas = etas,
                      lambdas = lambdas)

write_tsv(cones, 'cones.tsv')
