#!/usr/bin/env nextflow

params.out = '.'

params.network_permutations = 100
params.heat_permutations = 1000
beta = 0.4
params.fdr = 0.2

process make_network {

    input:
        path EDGELIST
        val BETA
        val NETWORK_PERMUTATIONS

    output:
        tuple path("ppin_ppr_${BETA}.h5"), path('permuted')

    """
    R --no-save <<code
library(tidyverse)
library(igraph)

net <- read_tsv("${EDGELIST}") %>%
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
        val CUTOFF

    output:
        path "${SCORES.getBaseName()}.json"

    """
    R --no-save <<code
library(tidyverse)

read_tsv('${SCORES}') %>%
    # sparsify scores
    mutate(padj = p.adjust(Pvalue, method = "fdr"),
           p = ifelse(padj < ${CUTOFF}, Pvalue, 1), 
           score = -log10(p)) %>%
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
--output_directory .

    mv consensus/subnetworks.tsv ${HEAT.getBaseName()}.subnetworks.tsv
    """

}

process process_output {

    publishDir params.out, mode: 'copy'
    tag { SUBNETWORKS.getBaseName() }

    input:
        path SCORES
        path SUBNETWORKS

    output:
        path "${SCORES.getBaseName()}.hotnet2.tsv"

    """
#!/usr/bin/env Rscript

library(tidyverse)

read_tsv('${SUBNETWORKS}', col_types = 'cc', skip = 1) %>%
    rename(gene = `#Core`) %>%
    select(gene) %>%
    mutate(cluster = ifelse(n() > 0, 1:n(), 0)) %>%
    separate_rows(gene, sep = ' ') %>%
    write_tsv('${SCORES.getBaseName()}.hotnet2.tsv')
    """

}

workflow hotnet2_nf {
    take:
        scores
        edgelist
        fdr
        network_permutations
        heat_permutations
    main:
        make_network(edgelist, beta, network_permutations)
        compute_heat(scores, fdr)
        hotnet2(make_network.out, compute_heat.out, beta, network_permutations, heat_permutations)
        process_output(scores, hotnet2.out)
    emit:
        process_output.out
}

workflow {
    main:
        hotnet2_nf(file(params.scores), file(params.edgelist), params.fdr, params.network_permutations, params.heat_permutations)
    emit:
        hotnet2_nf.out
}
