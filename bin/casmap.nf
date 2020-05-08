#!/usr/bin/env nextflow

params.out = '.'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

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

process casmap {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file RGWAS from rgwas

    output:
        file 'scored_interactions.casmap.tsv'

    """
    #!/usr/bin/env Rscript

    library(CASMAP)
    library(tidyverse)

    load('${RGWAS}')

    # make data files
    X <- as(gwas[['genotypes']], 'numeric')
    X[is.na(X)] = 0 # change?
    X[X >= 1] = 1 # dominant encoding
    write.table(t(X), 'genotypes.tsv', sep = ' ',
                col.names = FALSE, row.names = FALSE)

    y <- gwas[['fam']][['affected']] - 1
    write.table(y, 'phenotypes.tsv', sep = ' ',
                col.names = FALSE, row.names = FALSE)

    # run casmap
    cm <- CASMAP(mode="higherOrderEpistasis")
    cm\$readFiles(genotype_file   = 'genotypes.tsv',
                  phenotype_file  = 'phenotypes.tsv') 
                  # covariate_file  = 'covariates.tsv')

    cm\$execute()

    # prepare output
    results <- cm\$getSignificantInteractions()
    snps <- rownames(gwas[['map']])

    mutate(results,
           itemsets = lapply(itemsets, function(i) snps[i]) %>% lapply(paste, collapse = ',') %>% unlist) %>% 
        write_tsv(path = 'scored_interactions.casmap.tsv')
    """

}