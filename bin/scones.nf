#!/usr/bin/env nextflow

include { get_bfile } from './templates/utils.nf'

params.out = '.'

// gwas
bfile = get_bfile(params.bfile)

// SConES parameters
params.network = 'gs'
params.score = 'chi2'
params.criterion = 'stability'
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

    publishDir params.out, mode: 'copy'
    tag { RGWAS.getBaseName() }

    input:
        tuple path(RGWAS), path(RNET)
        val SCORE
        val CRITERION

    output:
        path "${RGWAS.getBaseName()}.scones.tsv"

    """
#!/usr/bin/env Rscript
library(igraph)
library(martini)
library(tidyverse)

load("${RGWAS}")
load("${RNET}")

# make a exploratory run to get the best parameters
params <- capture.output(
    scones.cv(gwas, net, score = "${SCORE}", criterion = "${CRITERION}"),
    type = "message")

cat(params)

params <- tail(params, n = 2) %>% 
  lapply(strsplit, ' =') %>% 
  unlist %>% .[c(F,T)] %>% 
  as.numeric() %>% 
  log10

# optimize the parameters
etas <- 10^seq(params[1] - 1.5, params[1] + 1.5, length.out = 10)
lambdas <- 10^seq(params[2] - 1.5, params[2] + 1.5, length.out = 10)

cones <- scones.cv(gwas, net,
                   score = "${SCORE}",
                   criterion = "${CRITERION}",
                   etas = etas,
                   lambdas = lambdas)

tibble(snp = names(V(cones))) %>%
    write_tsv('${RGWAS.getBaseName()}.scones.tsv')
    """

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

    """
#!/usr/bin/env Rscript
library(martini)
library(igraph)
library(tidyverse)

load("$RGWAS")
load("$RNET")

cones <- search_cones(gwas, 
                      net,
                      associationScore = "${SCORE}",
                      modelScore = "${CRITERION}",
search_cones(gwas, 
             net,
             associationScore = "${SCORE}",
             modelScore = "${CRITERION}",
             etas = c(${ETA}),
             lambdas = c(${LAMBDA})) %>%
    filter(selected) %>%
    write_tsv('${RGWAS.getBaseName()}.scones.tsv')
    """

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

workflow {
    main:
        scones_nf(bfile, file(params.edgelist), params.network, file(params.snp2gene), params.score, params.criterion, params.eta, params.lambda)
    emit:
        scones_nf.out
}
