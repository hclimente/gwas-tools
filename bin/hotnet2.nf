#!/usr/bin/env nextflow

params.out = "."

tab2 = file(params.tab2)
vegas = file(params.scores)
hotnet2_path = file(params.hotnet2_path)

network_permutations = 100
heat_permutations = 1000
beta = 0.4

process make_network {

    input:
        file TAB2 from tab2

    output:
        file 'node_index.tsv' into node_index
        file 'edge_list.tsv' into edge_list

    script:
    template 'io/tab2_2hotnet.R'

}

process vegas2hotnet {

    input:
        file VEGAS from vegas

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
--output_file heat.json
    """

}

process hotnet2 {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file HOTNET2 from hotnet2_path
        file HEAT from heat
        file NETWORK from h5
        file PERMS from permutations
        val BETA from beta

    output:
        file 'selected_genes.hotnet2.txt'

    """
    python2 ${HOTNET2}/HotNet2.py \
--network_files ${NETWORK} \
--permuted_network_path ${PERMS}/ppin_ppr_${BETA}_##NUM##.h5 \
--heat_files ${HEAT} \
--network_permutations ${network_permutations} \
--heat_permutations ${heat_permutations} \
--output_directory .

    cp consensus/subnetworks.tsv selected_genes.hotnet2.txt
    """

}
