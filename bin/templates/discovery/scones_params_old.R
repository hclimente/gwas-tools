#!/usr/bin/env Rscript
library(martini)
library(igraph)
library(tidyverse)

load("$RGWAS")
load("$RNET")

cones <- search_cones(gwas, net,
                      associationScore = "${SCORE}",
                      modelScore = "${CRITERION}",
                      etas = c(${ETA}),
                      lambdas = c(${LAMBDA}))

write_tsv(cones, 'cones.tsv')
