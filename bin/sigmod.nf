#!/usr/bin/env nextflow

params.out = '.'

// sigmod params
SIGMOD_PATH = file(params.sigmod)
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10

// annotation
VEGAS_OUT = file(params.vegas)
TAB2 = file(params.tab2)

process sigmod {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file VEGAS_OUT
        file TAB2
        file SIGMOD_PATH
        val LAMBDAMAX from params.lambdamax
        val NMAX from params.nmax
        val MAXJUMP from params.maxjump

    output:
        file 'selected_genes.sigmod.txt'

    script:
    template 'discovery/run_sigmod.R'
}
