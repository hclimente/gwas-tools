#!/usr/bin/env nextflow

params.out = '.'

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

    publishDir params.out, mode: 'copy'
    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path EDGELIST
        val LAMBDAMAX
        val NMAX
        val MAXJUMP
        path SIGMOD_PATH

    output:
        path "${SCORES.getBaseName()}.sigmod.txt"

    script:
    template 'discovery/sigmod.R'

}

workflow sigmod_nf {
    take:
        scores
        edgelist
        lambdamax
        nmax
        maxjump
    main:
        download_sigmod()
        sigmod(scores, edgelist, lambdamax, nmax, maxjump, download_sigmod.out)
    emit:
        sigmod.out
}

workflow {
    main:
        sigmod_nf(file(params.scores), file(params.edgelist), params.lambdamax, params.nmax, params.maxjump)
    emit:
        sigmod_nf.out
}
