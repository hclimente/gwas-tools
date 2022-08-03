#!/usr/bin/env nextflow

params.out = "."

tab2 = file(params.tab2)
vegas = file(params.scores)
hotnet2_path = file(params.hotnet2_path)

network_permutations = 100
heat_permutations = 1000
beta = 0.4
params.lfdr_cutoff = 0.05

process make_network {

    input:
        file TAB2 from tab2

    output:
        file 'node_index.tsv' into node_index
        file 'edge_list.tsv' into edge_list

    script:
    template 'io/tab2_2hotnet.R'

}

process sparse_scores {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
    file SCORES from vegas
    val CUTOFF from params.lfdr_cutoff

    output:
    file 'scored_genes.sparse.txt' into sparse_scores
    file 'lfdr_plot.pdf'

    """
#!/usr/bin/env Rscript

library(tidyverse)
library(twilight)
library(cowplot)

theme_set(theme_cowplot())

scores <- read_tsv('${SCORES}')

lfdr <- twilight(scores\$Pvalue, B=1000)
lfdr <- tibble(Gene = scores\$Gene[as.numeric(rownames(lfdr\$result))],
               vegas_p = scores\$Pvalue[as.numeric(rownames(lfdr\$result))],
               lfdr = lfdr\$result\$fdr)

ggplot(lfdr, aes(x = vegas_p, y = 1 - lfdr)) +
    geom_line() +
    geom_vline(xintercept = ${CUTOFF}, color = 'red') +
    labs(x = 'P-value', y = '1 - lFDR')
ggsave('lfdr_plot.pdf', width=7, height=6)

lfdr %>%
    mutate(Pvalue = ifelse(vegas_p < ${CUTOFF}, vegas_p, 1)) %>%
    write_tsv('scored_genes.sparse.txt')
    """

}

process vegas2hotnet {

    input:
        file VEGAS from sparse_scores 

    output:
        file 'scores.ht' into scores

    script:
    template 'io/vegas2hotnet.R'

}

process make_h5_network {

    input:
        file HOTNET2 from hotnet2_path
        file NODE_IDX from node_index
        file EDGE_LIST from edge_list
        val BETA from beta

    output:
        file "ppin_ppr_${BETA}.h5" into h5
        file 'permuted' into permutations

    """
    python2 ${HOTNET2}/makeNetworkFiles.py \
--edgelist_file ${EDGE_LIST} \
--gene_index_file ${NODE_IDX} \
--network_name ppin \
--prefix ppin \
--beta ${BETA} \
--cores -1 \
--num_permutations ${network_permutations} \
--output_dir .
    """

}

process make_heat_data {

    input:
        file HOTNET2 from hotnet2_path
        file SCORES from scores

    output:
        file 'heat.json' into heat

    """
    python2 ${HOTNET2}/makeHeatFile.py \
scores \
--heat_file ${SCORES} \
--output_file heat.json \
--name gwas
    """

}

process hotnet2 {

    input:
        file HOTNET2 from hotnet2_path
        file HEAT from heat
        file NETWORK from h5
        file PERMS from permutations
        val BETA from beta

    output:
        file 'consensus/subnetworks.tsv' into subnetworks

    """
    python2 ${HOTNET2}/HotNet2.py \
--network_files ${NETWORK} \
--permuted_network_path ${PERMS}/ppin_ppr_${BETA}_##NUM##.h5 \
--heat_files ${HEAT} \
--network_permutations ${network_permutations} \
--heat_permutations ${heat_permutations} \
--num_cores -1 \
--output_directory .
    """

}

process process_output {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file SUBNETWORKS from subnetworks

    output:
        file 'selected_genes.hotnet2.tsv'

    """
#!/usr/bin/env Rscript

library(tidyverse)

read_tsv('${SUBNETWORKS}', col_types = 'cc', comment = '#', col_names = F) %>%
    select(X1) %>%
    mutate(cluster = 1:n()) %>%
    separate_rows(X1, sep = ' ') %>%
    rename(gene = X1) %>%
    write_tsv('selected_genes.hotnet2.tsv')
    """

}
