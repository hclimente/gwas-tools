#!/usr/bin/env nextflow

params.out = '.'

SIGMOD_PATH = file(params.sigmod)

// annotation
VEGAS_OUT = file(params.vegas)
TAB2 = file(params.tab2)

process run_sigmod {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file VEGAS_OUT
        file TAB2
        file SIGMOD_PATH

    output:
        file 'selected_genes.sigmod.txt'

    script:
    template 'discovery/run_sigmod.R'
}