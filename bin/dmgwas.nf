#!/usr/bin/env nextflow

params.r = 0.1
params.d = 2

process dmgwas {

    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path EDGELIST
        val D
        val R

    output:
        file 'selected_genes.dmgwas.txt'

    """
    #!/usr/bin/env Rscript

    library(dmGWAS)
    library(tidyverse)

    scores <- read_tsv('${SCORES}') %>% 
        select(Gene, Pvalue) %>%
        mutate(Pvalue = ifelse(Pvalue == 1, 0.99999, Pvalue)) %>%
        as.data.frame
    net <- read_tsv('${EDGELIST}')

    modules <- dms(net, scores, expr1 = NULL, expr2 = NULL, r = ${R}, d = ${D})
    top <- simpleChoose(modules)

    tibble(gene = names(V(top\$subnetwork))) %>%
        write_tsv('selected_genes.dmgwas.txt')
    """

}

workflow dmgwas_nf {
    take:
        scores
        tab2
        d
        r
    main:
        dmgwas(scores, tab2, d, r)
    emit:
        dmgwas.out
}

workflow {
    main:
        dmgwas_nf(file(params.scores), file(params.tab2), params.d, params.r)
    emit:
        dmgwas_nf.out
}
