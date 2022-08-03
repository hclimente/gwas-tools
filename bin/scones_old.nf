#!/usr/bin/env nextflow

include { read_bfile; make_snp_network } from './scones.nf'
include { get_bfile } from './templates/utils.nf'

// gwas
bfile = get_bfile(params.bfile)

// SConES parameters
params.network = 'gs'
params.score = 'chi2'
params.criterion = 'consistency'
params.encoding = 'additive'
params.eta = null
params.lambda = null

process scones {

    input:
        file RGWAS
        file RNET
        val SCORE
        val CRITERION

    output:
        file 'cones.tsv'

    script:
    template 'discovery/scones_old.R'

}

process parametrized_scones {

    input:
        file RGWAS
        file RNET
        val SCORE
        val CRITERION
        val ETA
        val LAMBDA

    output:
        path 'cones.tsv'

    script:
    template 'discovery/scones_params_old.R'

}

workflow scones_old_nf {
    take:
        bfile
        tab2
        network
        snp2gene
        score
        criterion
        eta
        lambda
    main:
        snp2gene = (network == 'gm' | network == 'gi') ? file(snp2gene) : file('NO_SNP2GENE')
        tab2 = (network == 'gi') ? file(tab2) : file('NO_TAB2')

        read_bfile(bfile)
        make_snp_network(tab2, network, snp2gene, read_bfile.out)

        if (eta == null | lambda == null) {
            scones(read_bfile.out, make_snp_network.out, score, criterion)
            out = scones.out
        } else {
            parametrized_scones(read_bfile.out, make_snp_network.out, score, criterion, eta, lambda)
            out = parametrized_scones.out
        }
    emit:
        out
}

workflow {
    main:
        scones_old_nf(bfile, params.tab2, params.network, snp2gene, params.score, params.criterion, params.eta, params.lambda)
    emit:
        scones_old_nf.out
}
