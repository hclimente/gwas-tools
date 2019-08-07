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

process gpuepiscan {

        publishDir "$params.out", overwrite: true, mode: "copy"

        input:
          file RGWAS from rgwas

        output:
          file 'scored_interactions.episcan.tsv'

        script:
        """
        #!/usr/bin/env Rscript

        library(gpuEpiScan)
        library(tidyverse)
        library(bigmemory)

        load('${RGWAS}')
        
        # read dataset
        X <- as(gwas[['genotypes']], 'numeric') # %>% as.big.matrix
        y <- gwas[['fam']][['affected']] - 1

        gpuEpiScan(geno1 = X,
                   pheno = y,
                   phetype = "quantitative",
                   outfile = "scored_interactions.episcan", 
                   suffix = ".txt" 
                  )
        
        system2('sed', args = c("'s/ /\t/g'", 'scored_interactions.episcan.txt'),
                stdout = 'scored_interactions.episcan.tsv')
        """

}
