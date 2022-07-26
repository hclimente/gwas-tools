#!/usr/bin/env nextflow

params.out = '.'
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10
params.sigmod = null
docker_sigmod = '/gwas-tools/SigMod_v2'

// conditional SigMod input to handle Docker or Dockerless execution
if (params.sigmod != null) {
    SIGMOD_PATH = file(params.sigmod) 
} else {
    SIGMOD_PATH = docker_sigmod
}

// annotation
VEGAS_OUT = file(params.vegas)
TAB2 = file(params.tab2)

process sigmod {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file VEGAS_OUT
        file TAB2
        if (params.sigmod != null) {
           file SIGMOD_PATH
        } else {
            val SIGMOD_PATH
        }
        val LAMBDAMAX from params.lambdamax
        val NMAX from params.nmax
        val MAXJUMP from params.maxjump

    output:
        file 'selected_genes.sigmod.txt'

    script:
    template 'discovery/run_sigmod.R'
}
