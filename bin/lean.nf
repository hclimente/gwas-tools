#!/usr/bin/env nextflow

params.out = '.'

// annotation
VEGAS_OUT = file(params.vegas)
TAB2 = file(params.tab2)

process tab22igraph {

    input:
        file TAB2
    
    output:
        file 'net.RData' into RNET

    script:
    template 'io/tab22igraph.R'

}

process read_vegas {

    input:
        file VEGAS_OUT

    output:
        file 'scores.RData' into RSCORES

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    vegas_out <- read_tsv('$VEGAS_OUT')

    scores <- vegas_out[['Pvalue']]
    names(scores) <- vegas_out[['Gene']]

    save(scores, file = 'scores.RData')
    """

}

process lean {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file RSCORES
        file RNET

    output:
        file 'scored_genes.lean.txt'

    script:
    template 'discovery/run_leanr.R'
}
