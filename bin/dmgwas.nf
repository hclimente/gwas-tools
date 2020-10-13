#!/usr/bin/env nextflow

params.out = '.'
params.r = 0.1
params.d = 2

// annotation
VEGAS_OUT = file(params.vegas)
TAB2 =  file(params.tab2)

process dmgwas {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file VEGAS_OUT
        file TAB2
        val D from params.d
        val R from params.r

    output:
        file 'selected_genes.dmgwas.txt'

    """
    #!/usr/bin/env Rscript

    library(dmGWAS)
    library(tidyverse)

    vegas <- read_tsv('${VEGAS_OUT}') %>% 
        select(Gene, Pvalue) %>%
        mutate(Pvalue = ifelse(Pvalue == 1, 0.99999, Pvalue)) %>%
        as.data.frame
    ppi <- read_tsv('${TAB2}',
		    col_types = 'cccccccccccccccccccccccc') %>%
        rename(interactorA = `Official Symbol Interactor A`, 
               interactorB = `Official Symbol Interactor B`) %>%
        select(interactorA, interactorB)

    modules <- dms(ppi, vegas, expr1 = NULL, expr2 = NULL, r = ${R}, d = ${D})
    top <- simpleChoose(modules)

    tibble(gene = names(V(top\$subnetwork))) %>%
        write_tsv('selected_genes.dmgwas.txt')
    """

}
