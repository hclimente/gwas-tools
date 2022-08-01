#!/usr/bin/env nextflow

process lean {

    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path TAB2

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
        lean(scores, tab2)
    emit:
        lean.out
}

workflow {
    main:
        lean_nf(file(params.scores), file(params.tab2))
    emit:
        lean_nf.out
}
