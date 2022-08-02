#!/usr/bin/env nextflow

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")
bfile = tuple(bed, bim, fam)

// SConES parameters
params.network = 'gs'
params.score = 'chi2'
params.criterion = 'consistency'
params.encoding = 'additive'

// additional files
snp2gene = (params.network == 'gm' | params.network == 'gi') ? file(params.snp2gene) : file('NO_SNP2GENE')
tab2 = (params.network == 'gi') ? file(params.tab2) : file('NO_TAB2')

process bed2r {

    input:
        tuple path(BED), path(BIM), path(FAM)

    output:
        path 'gwas.RData'

    script:
    template 'io/bed2r.R'

}

process create_snp_network {

    input:
        file TAB2
        val NET
        file SNP2GENE
        file RGWAS

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
            tab2 <- read_tsv("${TAB2}") %>%
                rename(gene1 = `Official Symbol Interactor A`, gene2 = `Official Symbol Interactor B`) %>%
                select(gene1, gene2)
            net <- get_GI_network(gwas, snpMapping = snp2gene, ppi = tab2)
        }
    } else {
        stop("network type not recognized.")
    }

    save(net, file = "net.RData")
    """
}

process scones {

    input:
        file RGWAS
        file RNET
        val SCORE
        val CRITERION

    output:
        file 'cones.tsv'

    script:
    template 'discovery/scones.R'

}

workflow scones_nf {
    take:
        bfile
        tab2
        network_type
        snp2gene
        score
        criterion
    main:
        bed2r(bfile)
        make_snp_network(tab2, network_type, snp2gene, bed2r.out)
        scones(bed2r.out, make_snp_network.out, score, criterion)
    emit:
        scones.out
}

workflow {
    main:
        scones_nf(bfile, file(params.tab2), params.network, snp2gene, params.score, params.criterion)
    emit:
        scones_nf.out
}
