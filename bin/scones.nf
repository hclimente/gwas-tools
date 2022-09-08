#!/usr/bin/env nextflow

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

process make_gs_network {

    tag { RGWAS.getBaseName() }

    input:
        path RGWAS

    output:
        tuple path(RGWAS), path("gs.RData")

    """
    #!/usr/bin/env Rscript
    library(martini)

    load("${RGWAS}")

    net <- get_GS_network(gwas)
    save(net, file = "gs.RData")
    """

}

process make_gm_network {

    tag { RGWAS.getBaseName() }

    input:
        path RGWAS
        path SNP2GENE

    output:
        tuple path(RGWAS), path("gm.RData")

    """
    #!/usr/bin/env Rscript
    library(martini)
    library(tidyverse)

    load("${RGWAS}")

    net <- get_GM_network(gwas,
                          snpMapping = read_tsv("${SNP2GENE}"))
    save(net, file = "gm.RData")
    """

}

process make_gi_network {

    tag { RGWAS.getBaseName() }

    input:
        path RGWAS
        path SNP2GENE
        path EDGELIST

    output:
        tuple path(RGWAS), path("gi.RData")

    """
    #!/usr/bin/env Rscript
    library(martini)
    library(tidyverse)
    library(igraph)

    load("${RGWAS}")
    
    net <- get_GI_network(gwas,
                          snpMapping = read_tsv("${SNP2GENE}"),
                          ppi = read_tsv("${EDGELIST}"))
    save(net, file = "gi.RData")
    """

}

process scones {


    tag { RGWAS.getBaseName() }

    input:
        tuple path(RGWAS), path(RNET)
        val SCORE
        val CRITERION

    output:
        path "${RGWAS.getBaseName()}.scones.tsv"

    script:
    template 'discovery/scones.R'

}

process parametrized_scones {

    tag { RGWAS.getBaseName() }

    input:
        tuple path(RGWAS), path(RNET)
        val SCORE
        val CRITERION
        val ETA
        val LAMBDA

    output:
        path "${RGWAS.getBaseName()}.scones.tsv"

    script:
    template 'discovery/scones_params.R'

}

workflow scones_nf {
    take:
        bfile
        edgelist
        network
        snp2gene
        score
        criterion
        eta
        lambda
    main:
        read_bfile(bfile)

        switch(network) {
            case 'gs':
                make_gs_network(read_bfile.out)
                net = make_gs_network.out
                break
            case 'gm':
                net = make_gm_network(read_bfile.out, snp2gene)
                net = make_gm_network.out
                break
            case 'gi':
                net = make_gi_network(read_bfile.out, snp2gene, edgelist)
                net = make_gi_network.out
                break
        }


        if (eta == null | lambda == null) {
            scones(net, score, criterion)
            out = scones.out
        } else {
            parametrized_scones(net, score, criterion, eta, lambda)
            out = parametrized_scones.out
        }
    emit:
        out
}

workflow scones_new_nf {
    take:
        bfile
        edgelist
        network
        snp2gene
        score
        criterion
    main:
        snp2gene = (network == 'gm' | network == 'gi') ? snp2gene : null
        edgelist = (network == 'gi') ? edgelist : null

        read_bfile(bfile)
        make_snp_network(edgelist, network, snp2gene, read_bfile.out)
        scones(read_bfile.out, make_snp_network.out, score, criterion)
    emit:
        scones.out
}

workflow {
    main:
        scones_nf(bfile, file(params.edgelist), params.network, file(params.snp2gene), params.score, params.criterion, params.eta, params.lambda)
    emit:
        scones_nf.out
}
