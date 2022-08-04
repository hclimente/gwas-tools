#!/usr/bin/env nextflow

params.out = "."

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

process compute_grm {

    input:
    file BED from bed
    file BIM from bim
    file FAM from fam

    output:
    file 'plink.grm.gz' into grm

    """
plink -bfile ${BED.baseName} -make-grm-gz
"""
}

process lmm {

    input:
    file BED from bed
    file BIM from bim
    file FAM from fam
    file GRM from grm

    script:
    """
#!/usr/bin/env Rscript

library(GMMAT)
library(tidyverse)

load('${GRM}')

pheno <- read_tsv('${FAM}', col_names = FALSE) %>%
    rename(sex = X5, status = X6)

model0 <- glmmkin(disease ~ sex, data = pheno, kins = GRM,
                  id = "id", family = binomial(link = "logit"))

glmm.score(model0, infile = '${BED.baseName}', outfile = 'scored_snps.glmm.txt")

    """

}
