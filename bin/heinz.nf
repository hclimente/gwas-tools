#!/usr/bin/env nextflow

params.out = '.'
params.fdr = 0.1

process heinz {

    publishDir params.out, mode: 'copy'
    tag { SCORES.getBaseName() }

    input:
        file SCORES
        file EDGELIST
        val FDR

    output:
        path "${SCORES.getBaseName()}.heinz.txt"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(igraph)
    library(BioNet)

    vegas <- read_tsv('${SCORES}')
    net <- read_tsv("${EDGELIST}") %>%
        filter(gene1 %in% vegas\$Gene & gene2 %in% vegas\$Gene) %>%
        graph_from_data_frame(directed = FALSE)
    vegas <- filter(vegas, Gene %in% names(V(net)))

    # search subnetworks
    pvals <- vegas\$Pvalue
    names(pvals) <- vegas\$Gene
    fb <- fitBumModel(pvals, plot = FALSE)
    scores <- scoreNodes(net, fb, fdr = ${FDR})

    if (sum(scores > 0)) {
        selected <- runFastHeinz(net, scores)    
    	tibble(gene = names(V(selected))) %>% 
            write_tsv("${SCORES.getBaseName()}.heinz.txt")
    } else {
        tibble(gene = character()) %>% 
            write_tsv("${SCORES.getBaseName()}.heinz.txt")
    }
    """

}

workflow heinz_nf {
    take:
        scores
        edgelist
        fdr
    main:
        heinz(scores, edgelist, fdr)
    emit:
        heinz.out
}

workflow {
    main:
        heinz_nf(file(params.scores), file(params.edgelist), params.fdr)
    emit:
        heinz_nf.out
}
