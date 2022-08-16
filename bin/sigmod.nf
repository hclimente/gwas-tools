#!/usr/bin/env nextflow

// sigmod params
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10

process download_sigmod {

    output:
        path "SigMod_v2"

    """
    wget https://github.com/YuanlongLiu/SigMod/raw/20c561876d87a0faca632a6b93882fcffd719b17/SigMod_v2.zip && unzip SigMod_v2.zip
    """
}

process sigmod {

    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path EDGELIST
        val LAMBDAMAX
        val NMAX
        val MAXJUMP
        path SIGMOD_PATH

    output:
        path 'selected_genes.sigmod.txt'

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
        download_sigmod()
        sigmod(scores, tab2, lambdamax, nmax, maxjump, download_sigmod.out)
    emit:
        sigmod.out
}

workflow {
    main:
        sigmod_nf(file(params.scores), file(params.tab2), params.lambdamax, params.nmax, params.maxjump)
    emit:
        sigmod_nf.out
}
