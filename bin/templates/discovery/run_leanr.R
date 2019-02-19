#!/usr/bin/env Rscript
library(LEANR)
library(tidyverse)

load("${RSCORES}")
load("${RNET}")

scores <- scores[names(scores) %in% names(V(net))]

results <- run.lean(scores, net, n_reps = 10000, 
                    add.scored.genes = TRUE, verbose = TRUE)

write_tsv(as.data.frame(results\$restab), 'scored_genes.lean.txt')