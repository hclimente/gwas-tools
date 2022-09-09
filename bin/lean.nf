#!/usr/bin/env nextflow

process lean {

    publishDir params.out, mode: 'copy'
    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path EDGELIST

    output:
        path "${SCORES.getBaseName()}.lean.tsv"

    script:
    template 'discovery/leanr.R'
}

workflow lean_nf {
    take:
        scores
        edgelist
    main:
        lean(scores, edgelist)
    emit:
        lean.out
}

workflow {
    main:
        lean_nf(file(params.scores), file(params.edgelist))
    emit:
        lean_nf.out
}
