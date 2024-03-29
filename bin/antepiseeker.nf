#!/usr/bin/env nextflow

params.out = "."

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

n_snps = Channel
        .fromPath("${bim}")
        .splitText()
        .count()

process bed2r {

    input:
        file BED from bed
        file BIM from bim
        file FAM from fam

    output:
        file 'gwas.RData' into rdata

    script:
    template 'io/bed2r.R'

}

process r2aes {

    input:
        file rdata

    output:
        file "genotypes.aes" into aes_data

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(snpStats)
    load("$rdata")

    X <- as(gwas\$genotypes, "numeric")
    Y <- gwas\$fam\$affected - 1

    cbind(X,Y) %>% as.data.frame %>% write_csv("genotypes.aes")
    """

}

process antepiseeker {

    publishDir "$params.out", overwrite: true, mode: "copy"
    //validExitStatus 0,1

    input:
        val N from n_snps
        file AES_IN from aes_data

    output:
        file "scored_interactions.antepiseeker.txt"

    """
    cat << EOF >parameters.txt
    iAntCount        5000         // number of ants
    iItCountLarge    ${(N * 0.1).intValue()}   // number of iterations for the large haplotypes
    iItCountSmall    ${(N * 0.05).intValue()}  // number of iterations for the small haplotypes
    alpha            1            // weight given to pheromone deposited by ants
    iTopModel        1000         // number of top ranking haplotypes in the first stage
    iTopLoci         1000         // number of loci with top ranking pheromone in the first stage
    rou              0.05         // evaporation rate in Ant Colony Optimization
    phe              100          // initial pheromone level for each locus
    largehapsize     6            // size of the large haplotypes
    smallhapsize     3            // size of the small haplotypes
    iEpiModel        2            // number of SNPs in an epistatic interaction
    pvalue           0.01         // p value threshold (after Bonferroni correction)
    INPFILE          ${AES_IN}    // input file name for case-control genotype data
    OUTFILE          result.txt   // output file name for detected epistatic interactions
    EOF

    AntEpiSeeker

    echo -e "snp1\\tsnp2\\tchi2\\tp" >scored_interactions.antepiseeker.txt
    tail -n+3 result.txt | sed -E 's/[0-9]+\\(//g' | sed 's/)//g' | sed -E 's/[ \\t]+/\\t/g' >>scored_interactions.antepiseeker.txt
    """

}
