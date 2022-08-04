#!/usr/bin/env nextflow

params.out = '.'

// annotation
tab2 = file(params.tab2)

process rewire {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file TAB2 from tab2

    output:
        file 'rewired.tab2' 

    """
    #!/usr/bin/env Rscript

    #!/usr/bin/env Rscript
    library(igraph)
    library(tidyverse)

    read_tsv("${TAB2}") %>%
        rename(gene1 = `Official Symbol Interactor A`,
               gene2 = `Official Symbol Interactor B`) %>%
        select(gene1, gene2) %>%
        graph_from_data_frame(directed = FALSE) %>%
        rewire(with = keeping_degseq(niter = 300000)) %>%
        as_edgelist %>%
        as.data.frame %>%
        rename(`Official Symbol Interactor A` = V1, `Official Symbol Interactor B` = V2) %>%
        write_tsv('rewired.tab2')
    """

}
