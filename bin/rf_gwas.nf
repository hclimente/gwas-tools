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
