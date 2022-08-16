#!/usr/bin/env nextflow

include { get_bfile } from './templates/utils.nf'

// gwas
bfile = get_bfile(params.bfile)

// SConES parameters
params.snp2gene = null
params.tab2 = null

params.network = 'gs'
params.score = 'chi2'
params.criterion = 'consistency'
params.encoding = 'additive'
params.eta = null
params.lambda = null

process read_bfile {

    tag { BED.getBaseName() }
    afterScript "mv gwas.RData ${BED.getBaseName()}.Rdata"

    input:
        tuple path(BED), path(BIM), path(FAM)

    output:
        path "${BED.getBaseName()}.Rdata"

    script:
    template 'io/bed2r.R'

}

process make_snp_network {

    tag { RGWAS.getBaseName() }

    input:
        path EDGELIST
        val NET
        path SNP2GENE
        path RGWAS

    output:
        path 'net.RData'

    """
    #!/usr/bin/env Rscript
    library(martini)
    library(tidyverse)
    library(igraph)

    load("${RGWAS}")
    netType <- "${NET}"

    if (netType == "gs") {
        net <- get_GS_network(gwas)
    } else if (netType %in% c('gm', 'gi')) {
        snp2gene <- read_tsv("${SNP2GENE}")
        if (netType == "gm") {
            net <- get_GM_network(gwas, snpMapping = snp2gene)
        } else if (netType == "gi") {
            net <- get_GI_network(gwas, snpMapping = snp2gene, ppi = read_tsv("${EDGELIST}"))
        }
    } else {
        stop("network type not recognized.")
    }

    save(net, file = "net.RData")
    """
}

process scones {

    tag { RGWAS.getBaseName() }

    input:
        path RGWAS
        path RNET
        val SCORE
        val CRITERION

    output:
        path "${RGWAS.getBaseName()}.scones.tsv"

    script:
    template 'discovery/scones.R'

}

process scones_old {

    tag { RGWAS.getBaseName() }

    input:
        path RGWAS
        path RNET
        val SCORE
        val CRITERION

    output:
        path "${RGWAS.getBaseName()}.scones.tsv"

    script:
    template 'discovery/scones_old.R'

}

process parametrized_scones_old {

    tag { RGWAS.getBaseName() }

    input:
        path RGWAS
        path RNET
        val SCORE
        val CRITERION
        val ETA
        val LAMBDA

    output:
        path "${RGWAS.getBaseName()}.scones.tsv"

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
        snp2gene = (network == 'gm' | network == 'gi') ? snp2gene : null
        tab2 = (network == 'gi') ? tab2 : null

        read_bfile(bfile)
        make_snp_network(tab2, network, snp2gene, read_bfile.out)

        if (eta == null | lambda == null) {
            scones_old(read_bfile.out, make_snp_network.out, score, criterion)
            out = scones_old.out
        } else {
            parametrized_scones_old(read_bfile.out, make_snp_network.out, score, criterion, eta, lambda)
            out = parametrized_scones_old.out
        }
    emit:
        out
}

workflow scones_nf {
    take:
        bfile
        tab2
        network
        snp2gene
        score
        criterion
    main:
        snp2gene = (network == 'gm' | network == 'gi') ? snp2gene : null
        tab2 = (network == 'gi') ? tab2 : null

        read_bfile(bfile)
        make_snp_network(tab2, network, snp2gene, read_bfile.out)
        scones(read_bfile.out, make_snp_network.out, score, criterion)
    emit:
        scones.out
}

workflow {
    main:
        scones_old_nf(bfile, params.tab2, params.network, params.snp2gene, params.score, params.criterion, params.eta, params.lambda)
    emit:
        scones_old_nf.out
}
