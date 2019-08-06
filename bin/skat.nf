#!/usr/bin/env nextflow

params.out = "."
params.phenotype = 'discrete'

phenotype = (params.phenotype == 'discrete')? 'D' : 'C'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

// snp2gene
snp2gene = file(params.snp2gene)

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

process get_null {

    input:
    file RGWAS from rgwas
    val PHENO from phenotype

    output:
    file 'null.RData' into dnull

    """
    #!/usr/bin/env Rscript 

    library(SKAT)

    load('${RGWAS}')

    y <- gwas\$fam\$affected

    if ("${PHENO}" == 'D')
        y <- y - 1

    dnull <- SKAT_Null_Model(y ~ 1, out_type="${PHENO}")

    save(dnull, file = 'null.RData')
    """

}

process skat {

    publishDir "$params.out", overwrite: true, mode: "move"

    input:
    file RGWAS from rgwas
    file DNULL from dnull
    file SNP2GENE from snp2gene

    output:
    file 'scored_genes.skat.tsv'

    """
    #!/usr/bin/env Rscript 

    library(SKAT)
    library(tidyverse)
    library(snpStats)

    load('${RGWAS}')
    load('${DNULL}')

    snp2gene <- read_tsv('$SNP2GENE', col_types = 'cc')
    genes <- snp2gene\$gene %>% unique

    pvals <- lapply(genes, function(g){

        print(g)
        snps <- snp2gene\$snp[ snp2gene\$gene == g]
        Z <- as(gwas\$genotypes[,snps], 'numeric')
        Z <- abs(Z - 2)
        SKAT(Z, dnull)\$p.value

    }) %>% do.call(c, .)

    tibble(gene = genes, p = pvals) %>%
        write_tsv('scored_genes.skat.tsv')
    """

}
