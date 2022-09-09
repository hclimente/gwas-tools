#!/usr/bin/env Rscript
library(igraph)
library(LEANR)
library(tidyverse)

# read edgelist
net <- read_tsv("${EDGELIST}") %>%
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
    mutate(gene = rownames(results\$restab)) %>%
    select(gene, everything()) %>%
    write_tsv('${SCORES.getBaseName()}.lean.tsv')
