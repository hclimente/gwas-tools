#!/usr/bin/env Rscript
# TAB2
# node_index.tsv, edge_list.tsv
library(tidyverse)
library(igraph)

net <- read_tsv("${TAB2}") %>%
  rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
  select(gene1, gene2) %>%
  graph_from_data_frame(directed = FALSE)

as_edgelist(net, names = FALSE) %>% as.data.frame %>% write_tsv('edge_list.tsv', col_names = FALSE)
tibble(id = V(net), names = names(V(net))) %>% write_tsv('node_index.tsv', col_names = FALSE)
