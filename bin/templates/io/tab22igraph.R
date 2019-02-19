#!/usr/bin/env Rscript
library(igraph)
library(tidyverse)

net <- read_tsv("${TAB2}") %>%
    rename(gene1 = `Official Symbol Interactor A`, 
           gene2 = `Official Symbol Interactor B`) %>%
    select(gene1, gene2) %>%
    graph_from_data_frame(directed = FALSE)
save(net, file = 'net.RData')