#!/usr/bin/env Rscript

require(igraph)
library(tidyverse)

scripts <- list.files('${SIGMOD_PATH}/R', pattern='*.R\$', full.names=TRUE, ignore.case=TRUE)
sapply(scripts, source, .GlobalEnv)

# read network
net <- read_tsv("${EDGELIST}") %>%
    filter(gene1 != gene2) %>%
    as.data.frame

# read vegas output
scores <- read_tsv('${SCORES}') %>% 
    rename(gene = Gene, p = Pvalue) %>%
    select(gene, p) %>%
    as.data.frame

# check weight_index = NULL
scored_net <- construct_scored_net(net, interaction_indices = c(1,2), gene_ps = scores)
res_info <- SigMod_bisection(net = scored_net, lambda_max = ${LAMBDAMAX}, nmax = ${NMAX}, maxjump = ${MAXJUMP})

save(scored_net, res_info, file = 'sigmod.RData')

data.frame(gene = names(V(res_info\$opt_module[[1]]))) %>%
    write_tsv("${SCORES.getBaseName()}.sigmod.txt")
