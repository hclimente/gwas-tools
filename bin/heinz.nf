#!/usr/bin/env nextflow

params.out = '.'
params.fdr = 0.1
nextflow.enable.dsl=2

process heinz {

    // publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file VEGAS
        file TAB2
        val FDR

    output:
        file 'selected_genes.heinz.txt' 

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(igraph)
    library(BioNet)

    vegas <- read_tsv('${VEGAS}')
    net <- read_tsv("${TAB2}") %>%
        rename(gene1 = `Official Symbol Interactor A`, 
               gene2 = `Official Symbol Interactor B`) %>%
        filter(gene1 %in% vegas\$Gene & gene2 %in% vegas\$Gene) %>%
        select(gene1, gene2) %>%
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
            write_tsv('selected_genes.heinz.txt')
    } else {
        write_tsv(tibble(gene = character()), 'selected_genes.heinz.txt')
    }
    """

}

workflow heinz_workflow {
    take:
        scores
        tab2
        fdr
    main:
        heinz(scores, tab2, fdr)
    emit:
        heinz.out
}

workflow {
    main:
        heinz_workflow(file(params.scores), file(params.tab2), params.fdr)
    emit:
        heinz_workflow.out
}