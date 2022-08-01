#!/usr/bin/env nextflow

process tab22igraph {

    input:
        path TAB2
    
    output:
        path 'net.RData'

    script:
    template 'io/tab22igraph.R'

}

process read_vegas {

    tag { VEGAS_OUT.getBaseName() }

    input:
        path VEGAS_OUT

    output:
        path "${VEGAS_OUT.getBaseName()}.RData"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    vegas_out <- read_tsv('$VEGAS_OUT')

    scores <- vegas_out[['Pvalue']]
    names(scores) <- vegas_out[['Gene']]

    save(scores, file = "${VEGAS_OUT.getBaseName()}.RData")
    """

}

process lean {

    tag { RSCORES.getBaseName() }

    input:
        path RSCORES
        path RNET

    output:
        path 'scored_genes.lean.txt'

    script:
    template 'discovery/leanr.R'
}

workflow lean_nf {
    take:
        scores
        tab2
    main:
        read_vegas(scores)
        tab22igraph(tab2)
        lean(read_vegas.out, tab22igraph.out)
    emit:
        lean.out
}

workflow {
    main:
        lean_nf(file(params.scores), file(params.tab2))
    emit:
        lean_nf.out
}
