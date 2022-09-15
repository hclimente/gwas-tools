#!/usr/bin/env nextflow

params.out = '.'
params.r = 0.1
params.d = 2

process dmgwas {

    publishDir params.out, mode: 'copy'
    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path EDGELIST
        val D
        val R

    output:
        path "${SCORES.getBaseName()}.dmgwas.txt"

    """
    #!/usr/bin/env Rscript

    library(dmGWAS)
    library(readr)
    library(dplyr)

    scores <- read_tsv('${SCORES}') %>% 
        select(Gene, Pvalue) %>%
        mutate(Pvalue = ifelse(Pvalue == 1, 0.99999, Pvalue)) %>%
        as.data.frame
    net <- read_tsv('${EDGELIST}')

    modules <- dms(net, scores, expr1 = NULL, expr2 = NULL, r = ${R}, d = ${D})
    top <- simpleChoose(modules)

    tibble(gene = names(V(top\$subnetwork))) %>%
        write_tsv("${SCORES.getBaseName()}.dmgwas.txt")
    """

}

workflow dmgwas_nf {
    take:
        scores
        edgelist
        d
        r
    main:
        dmgwas(scores, edgelist, d, r)
    emit:
        dmgwas.out
}

workflow {
    main:
        dmgwas_nf(file(params.scores), file(params.edgelist), params.d, params.r)
    emit:
        dmgwas_nf.out
}
