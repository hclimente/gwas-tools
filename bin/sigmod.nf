#!/usr/bin/env nextflow

// sigmod params
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10

process sigmod {

    input:
        file SCORES
        file TAB2
        val LAMBDAMAX
        val NMAX
        val MAXJUMP

    output:
        file 'selected_genes.sigmod.txt'

    script:
    template 'discovery/sigmod.R'

}


workflow sigmod_nf {
    take:
        scores
        tab2
        lambdamax
        nmax
        maxjump
    main:
        sigmod(scores, tab2, lambdamax, nmax, maxjump)
    emit:
        sigmod.out
}

workflow {
    main:
        sigmod_nf(file(params.scores), file(params.tab2), params.lambdamax, params.nmax, params.maxjump)
    emit:
        sigmod_nf.out
}
