#!/usr/bin/env Rscript
# VEGAS
# scores.ht

library(tidyverse)

read_tsv('${VEGAS}') %>%
  mutate(score = -log10(Pvalue)) %>%
  select(Gene, score) %>%
  write_tsv('scores.ht', col_names = FALSE)
