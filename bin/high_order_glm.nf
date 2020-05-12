#!/usr/bin/env nextflow

params.out = '.'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

// tested snps
snp_set = file("${params.snp_set}")

process bed2r {

    input:
        file BED from bed
        file BIM from bim
        file FAM from fam

    output:
        file 'gwas.RData' into rgwas

    script:
    template 'io/bed2r.R'

}

process high_order_glm {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file RGWAS from rgwas
        file SNPS from snp_set

    output:
        file 'scored_interactions.high_order_glm.tsv'

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    snps <- read_tsv('${SNPS}')\$snp

    load('${RGWAS}')

    X <- as(gwas[['genotypes']], 'numeric')
    X[is.na(X)] = 0 # TODO change?
    X[X == 0] = 'AA'
    X[X == 1] = 'Aa'
    X[X == 2] = 'aa'

    y <- gwas[['fam']][['affected']] - 1

    intx <- model.matrix(~.^100, data = as.data.frame(X[, snps])) %>%
        as_tibble %>%
        mutate(y = y)

    fit <- glm(y ~ ., data = intx)
    pvals <- coef(summary(fit))[,'Pr(>|t|)']

    tibble(snps = names(pvals),
           p_value = pvals) %>%
        write_tsv('scored_interactions.high_order_glm.tsv')
    """

}
