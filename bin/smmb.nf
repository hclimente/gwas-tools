#!/usr/bin/env nextflow

params.out = "."

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

n_snps = Channel
        .fromPath("${BIM}")
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

process r2smmb {

    input:
        file rdata

    output:
        file "smmb_genotypes.txt" into smmb_genotypes
        file "smmb_phenotypes.txt" into smmb_phenotypes

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(snpStats)
    load("$rdata")

    y <- gwas\$fam\$affected - 1
    tibble(Class = y) %>% write_tsv('smmb_phenotypes.txt')

    X <- as(gwas\$genotypes, "numeric")
    as.tibble(X) %>% write_csv("smmb_genotypes.txt")
    """

}

process smmb {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file SMMB_GENOTYPES from smmb_genotypes
        file SMMB_PHENOTYPES from smmb_phenotypes

    output:
        file "scored_interactions.smmb.txt"

    """
    cat << EOF >PARAMETERS_SMMB.txt
    # number of header lines in the genotype file
    header 1

    # separation between fields in the genotype file
    separator ,

    # number of consecutive runs of SMMB (should be kept to this value)
    n_smmb_runs 1

    # alpha type I error rate
    alpha 0.05

    # precision used in adaptative permutations
    precision 0.03

    # size of the largest subset, K (draw in Algorithm 1)
    subset_size_large 10

    # size of the smallest subset, k (draw in Algorithm 2)
    subset_size_small 3

    # Number of Markov blankets learned to make the consensus
    n_mbs 100

    # Maximal number of iterations if there is no convergence (top level procedure of SMMB algorithm)
    n_trials_to_learn_mbs 100

    # Maximal number of iterations allowed to learn one Markov blanket, for one given large subset K, to escape the issue of non-modified MB.
    # If the sampling of k variables among K variables does not allow to identify a candidate subset s which is actually added to the growing MB, then the MB would not be modified and learnMB (inner level procedure of SMMB algorithm) would stop.
    # To palliate this problem, one coerces the exploration of the subset of K variables through n_trials_to_learn_1_mb.
    n_trials_to_learn_1_mb 30
    EOF

    SMMB ${SMMB_GENPTYPES} ${SMMP_PHENOTYPES}
    """

}