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
        file 'gwas.RData' into rgwas_train, rgwas_test

    script:
    template 'io/bed2r.R'

}

process random_forest {

    input:
        file RGWAS from rgwas_train

    output:
        file 'forest.RData' into forest
        file 'snp_pairs.tsv' into pairs

    """
    #!/usr/bin/env Rscript

    library(ranger)
    library(tidyverse)

    load("${RGWAS}")

    X <- as(gwas\$genotypes, "numeric")
    X[is.na(X)] <- 0 # safeguard against missing genotypes
    Y <- gwas\$fam\$affected

    gwas <- as_tibble(X) %>%
        mutate(phenotype = Y)

    rf <- ranger(dependent.variable.name = 'phenotype', data = gwas, importance = 'impurity')
    save(rf, file = 'forest.RData')

    snps <- lapply(seq(1, 500), function(i) {
        nodes <- treeInfo(rf,i)[['splitvarName']]
        names(nodes) <- treeInfo(rf,i)[['nodeID']]
        bind_rows(
                  treeInfo(rf,i)[,c('nodeID','leftChild')] %>% rename(otherNodeID = leftChild),
                  treeInfo(rf,i)[,c('nodeID','rightChild')] %>% rename(otherNodeID = rightChild)
                 ) %>%
            mutate(snp1 = nodes[nodeID + 1], snp2 = nodes[otherNodeID + 1]) %>%
            filter(!is.na(snp2)) %>%
            select(snp1, snp2)
    }) %>%
        do.call(rbind, .) %>% 
        unique %>%
        filter(snp1 != snp2)
    imp <- importance(rf)
    top_snps <- tail(sort(imp), n = floor(length(imp) * 0.2)) %>% names
    snps %>%
         filter(snp1 %in% top_snps | snp2 %in% top_snps) %>%
        write_tsv('snp_pairs.tsv')
    """

}

pairs. splitText( keepHeader: true, by: 10000, file: true ). set{ pairs_split }

process permute {

    input:
        each I from 1..10
        file RGWAS from rgwas_test
        file RF from forest
        file PAIRS from pairs_split

    output:
        file 'mse' into mse

    """
    #!/usr/bin/env Rscript

    library(ranger)
    library(tidyverse)

    load("${RGWAS}")
    load("${RF}")

    X <- as(gwas\$genotypes, "numeric")
    X[is.na(X)] <- 0 # safeguard against missing genotypes
    X <- as.tibble(X)
    Y <- gwas\$fam\$affected

    pairs <- read_tsv('${PAIRS}')

    lapply(seq(nrow(pairs)), function(i){

        snp1 <- as.character(pairs[i,1])
        snp2 <- as.character(pairs[i,2])

        # preserve main effects
        gwas <- mutate(X, phenotype = Y) %>%
            split(., .\$phenotype) %>%
            lapply(function(df){
                df[,snp1] <- df[sample(1:nrow(df)), snp1]
                df[,snp2] <- df[sample(1:nrow(df)), snp2]
                return(df)
            }) %>%
            bind_rows

        pred <- predict(rf, data = gwas)
        mse <- sum((gwas\$phenotype - pred[['predictions']])^2)/nrow(gwas)
        
        w_main <- tibble(snp1 = snp1, snp2 = snp2, i = ${I}, mse = mse, preserved = 'main')

        # preserve main effect and interaction
        gwas <- mutate(X, phenotype = Y) %>%
            split(., .\$phenotype) %>%
            lapply(function(df){
                perm <- sample(1:nrow(df))
                df[,snp1] <- df[perm, snp1]
                df[,snp2] <- df[perm, snp2]
                return(df)
            }) %>%
            bind_rows

        pred <- predict(rf, data = gwas)
        mse <- sum((gwas\$phenotype - pred[['predictions']])^2)/nrow(gwas)

        w_both <- tibble(snp1 = snp1, snp2 = snp2, i = ${I}, mse = mse, preserved = 'both')

        bind_rows(w_main, w_both)

    }) %>%
        bind_rows %>%
        write_tsv('mse', col_names = FALSE)
    """

}

process join_mse {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file "mse*" from mse. collect()

    output:
        file "mse.tsv"

    """
    echo 'snp1\tsnp2\ti\tmse\tpreserved' >mse.tsv
    cat mse* | sort >>mse.tsv
    """

}
