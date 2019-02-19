#!/usr/bin/env Rscript
library(LEANR)
library(tidyverse)

load("${RSCORES}")
load("${RNET}")

scores <- scores[names(scores) %in% names(V(net))]

results <- run.lean(scores, net, n_reps = 10000, 
                    add.scored.genes = TRUE, verbose = TRUE)

as_tibble(results\$restab) %>%
    mutate(Gene = rownames(results\$restab)) %>%
    select(Gene, everything()) %>%
    write_tsv('scored_genes.lean.txt')