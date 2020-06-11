#!/usr/bin/env nextflow

params.out = '.'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")
covars = file("${params.covars}")

// tested snps
interactions = file("${params.interactions}")

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

process make_sets {

    input:
        file INTERACTIONS from interactions
        val I from 2..interactions.countLines()

    output:
        file 'snp_set' into snp_sets

    """
    echo snp >snp_set
    sed '${I}q;d' ${INTERACTIONS} | sed 's/_/\\n/g' >>snp_set
    """

}

process high_order_glm {

    input:
        file RGWAS from rgwas
        file SNPS from snp_sets
        file COVARS from covars

    output:
        file 'scores.tsv' into scores

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    snps <- read_tsv('${SNPS}')\$snp
    covars <- read_tsv('${COVARS}')

    load('${RGWAS}')

    X <- as(gwas[['genotypes']], 'numeric')
    X[is.na(X)] = 0 # TODO change?
    X[X == 0] = 'AA'
    X[X == 1] = 'Aa'
    X[X == 2] = 'aa'

    y <- gwas[['fam']][['affected']] - 1

    rm(gwas)

    intx <- model.matrix(~.^100,data = as.data.frame(X[, snps])) %>%
        as_tibble %>%
        mutate(y = y) %>%
        bind_cols(covars)

    rm(X,y)

    fit <- glm(y ~ ., data = intx)
    pvals <- coef(summary(fit))[,'Pr(>|t|)']

    tibble(snp_set = paste(snps, collapse = '_'),
           test = names(pvals),
           p_value = pvals) %>%
        write_tsv('scores.tsv', col_names = FALSE)
    """

}

process join_results {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file "*" from scores.collect()

    output:
        file 'scored_interactions.high_order_glm.tsv'

    """
    echo "snp_set\ttest\tpvalue" >scored_interactions.high_order_glm.tsv 
    cat * >>scored_interactions.high_order_glm.tsv
    """

}