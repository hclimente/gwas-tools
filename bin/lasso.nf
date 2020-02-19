#!/usr/bin/env nextflow

params.out = '.'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

process lasso {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file BED from bed
        file FAM from fam
        file BIM from bim

    output:
        file 'scored_snps.lasso.tsv'

    """
    #!/usr/bin/env Rscript
    library(biglasso)
    library(snpStats)
    library(tidyverse)

    # read dataset
    gwas <- read.plink("${BED}", "${BIM}", "${FAM}")
    X <- as(gwas[['genotypes']], 'numeric')
    ## remove NAs
    X[is.na(X)] <- 0	
    X <- as.big.matrix(X) 
    y <- gwas[['fam']][['affected']] - 1
    names(y) <- gwas[['fam']][['pedigree']] %>% as.character

    rm(gwas)

    # train and evaluate classifier
    cvfit <- cv.biglasso(X, y, penalty = 'lasso', family = "binomial")

    tibble(snp = rownames(cvfit\$fit\$beta), 
	   beta = cvfit\$fit\$beta[,cvfit\$lambda == cvfit\$lambda.min][-1]) %>%
        write_tsv('scored_snps.lasso.tsv')
    """

}
