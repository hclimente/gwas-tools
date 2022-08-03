#!/usr/bin/env Rscript
library(martini)
library(igraph)
library(tidyverse)

load("$RGWAS")
load("$RNET")

subnet <- scones.cv(gwas, net,
                      score = "${SCORE}",
                      criterion = "${CRITERION}",
                      etas = c(${ETA}),
                      lambdas = c(${LAMBDA}))
                      
cones <- martini:::get_snp_modules(gwas, subnet)
cones[['c']] <- martini:::snp_test(gwas, covars, "${SCORE}")
cones[['selected']] <- cones[['snp']] %in% names(V(subnet))

write_tsv(cones, 'cones.tsv')
