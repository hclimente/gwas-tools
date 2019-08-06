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


y <- gwas['fam']['affected']
dnull <- SKAT_Null_Model(y ~ ., out_type="${PHENO}")

save(dnull, 'null.RData')
"""

}

process skat {

    input:
    file RGWAS from rgwas
    file DNULL from dnull
    file SNP2GENE from snp2gene

    output:
    file 'scored_genes.skat.tsv' into gene_scores

    """
#!/usr/bin/env Rscript 

library(SKAT)
library(tidyverse)

load('${RGWAS}')
load('${DNULL}')

snp2gene <- read_tsv('$SNP2GENE')
genes <- snp2gene\$gene %>% unique

pvals <- lapply(genes, function(g){

snps <- snp2gene\$snp[ snp2gene\$gene == g]
Z <- as(gwas['genotypes'], 'numeric')

SKAT(Z, dnull)\$p.value

}) %>% do.call(c, .)

tibble(gene = genes, p = pvals) %>%
write_tsv('scored_genes.skat.tsv')
    """

}
