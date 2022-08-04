#!/usr/bin/env nextflow

params.network_permutations = 100
params.heat_permutations = 1000
beta = 0.4
params.lfdr_cutoff = 0.05

process sparse_scores {

    input:
        path SCORES
        val CUTOFF

    output:
        path 'scored_genes.sparse.txt'

    """
#!/usr/bin/env Rscript

library(tidyverse)
library(twilight)

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

process make_network {

    input:
        path TAB2
        val BETA
        val NETWORK_PERMUTATIONS

    output:
        tuple path("ppin_ppr_${BETA}.h5"), path('permuted')

    """
    R --no-save <<code
library(tidyverse)
library(igraph)

net <- read_tsv("${TAB2}") %>%
  select(\\`Official Symbol Interactor A\\`, \\`Official Symbol Interactor B\\`) %>%
  graph_from_data_frame(directed = FALSE)

as_edgelist(net, names = FALSE) %>% 
    as.data.frame %>% 
    write_tsv('edge_list.tsv', col_names = FALSE)
tibble(id = V(net), names = names(V(net))) %>%
    write_tsv('node_index.tsv', col_names = FALSE)
code

    makeNetworkFiles.py \
--edgelist_file edge_list.tsv \
--gene_index_file node_index.tsv \
--network_name ppin \
--prefix ppin \
--beta ${BETA} \
--cores -1 \
--num_permutations ${NETWORK_PERMUTATIONS} \
--output_dir .
    """

}

process compute_heat {

    tag { SCORES.getBaseName() }

    input:
        path SCORES

    output:
        path "${SCORES.getBaseName()}.json"

    """
    R --no-save <<code
library(tidyverse)

read_tsv('${SCORES}') %>%
  mutate(score = -log10(Pvalue)) %>%
  select(Gene, score) %>%
  write_tsv('scores.ht', col_names = FALSE)
code

    makeHeatFile.py \
scores \
--heat_file scores.ht \
--output_file ${SCORES.getBaseName()}.json \
--name gwas
    """

}

process hotnet2 {

    tag { HEAT.getBaseName() }
    afterScript "mv consensus/subnetworks.tsv ${HEAT.getBaseName()}.subnetworks.tsv"

    input:
        tuple path(NETWORK), path(PERMS) 
        path HEAT
        val BETA
        val NETWORK_PERMUTATIONS
        val HEAT_PERMUTATIONS

    output:
        path "${HEAT.getBaseName()}.subnetworks.tsv"

    """
    HotNet2.py \
--network_files ${NETWORK} \
--permuted_network_path ${PERMS}/ppin_ppr_${BETA}_##NUM##.h5 \
--heat_files ${HEAT} \
--network_permutations ${NETWORK_PERMUTATIONS} \
--heat_permutations ${HEAT_PERMUTATIONS} \
--num_cores -1 \
--output_directory . \
--deltas 10
# TODO rm previous line
    """

}

process process_output {

    tag { SUBNETWORKS.getBaseName() }

    input:
        path SUBNETWORKS

    output:
        path 'selected_genes.hotnet2.tsv'

    """
#!/usr/bin/env Rscript

library(tidyverse)

read_tsv('${SUBNETWORKS}', col_types = 'cc', skip = 1) %>%
    rename(gene = `#Core`) %>%
    select(gene) %>%
    mutate(cluster = ifelse(n() > 0, 1:n(), 0)) %>%
    separate_rows(gene, sep = ' ') %>%
    write_tsv('selected_genes.hotnet2.tsv')
    """

}

workflow hotnet2_nf {
    take:
        scores
        tab2
        lfdr_cutoff
        network_permutations
        heat_permutations
    main:
        /* sparse_scores(scores, lfdr_cutoff) */
        make_network(tab2, beta, network_permutations)
        compute_heat(scores)
        hotnet2(make_network.out, compute_heat.out, beta, network_permutations, heat_permutations)
        process_output(hotnet2.out)
    emit:
        process_output.out
}

workflow {
    main:
        hotnet2_nf(file(params.scores), file(params.tab2), params.lfdr_cutoff, params.network_permutations, params.heat_permutations)
    emit:
        hotnet2_nf.out
}
