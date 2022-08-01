#!/usr/bin/env Rscript
library(igraph)
library(LEANR)
library(tidyverse)

# read tab2 file
net <- read_tsv("${TAB2}") %>%
    select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
    graph_from_data_frame(directed = FALSE)

# read gene scores
gene_scores <- read_tsv('$SCORES') %>%
    filter(Gene %in% names(V(net)))

scores <- gene_scores[['Pvalue']]
names(scores) <- gene_scores[['Gene']]

# run lean
results <- run.lean(scores, net, n_reps = 10000, 
                    add.scored.genes = TRUE, verbose = TRUE)

as_tibble(results\$restab) %>%
    mutate(Gene = rownames(results\$restab)) %>%
    select(Gene, everything()) %>%
    write_tsv('scored_genes.lean.txt')
