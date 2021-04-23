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

process bed2ped {

    input:
        file BED from bed
        file BIM from bim
        file FAM from fam

    output:
        file 'out.ped' into ped

    script:
    template 'io/bed2ped.R'

}

process pMDR {

    publishDir "$params.out", overwrite: true, mode: "copy"
    //validExitStatus 0,1

    input:
        val N from n_snps
        file PED from ped

    output:
        file "scored_interactions.pmdr.txt"

    """
    cat << EOF >parameters.txt
    ########################### Input ###########################
    # The dataset to be used
    INPUTFILE                          ${PED}
    # Input format. Merlin has 6 columns in ped file and an optional .dat file
    MERLIN_FORMAT                      NO
    # Specify the location of the optional dat file (contains the labels for each locus MERLIN_DAT                         
    # Set the value associated with affected status
    AFFECTED_VALUE                     2
    # Set the value associated with the unaffected status
    # Individuals with neither Affected nor Unaffected will be ignored by the analysis
    UNAFFECTED_VALUE                   1
    # Ignore 0 or more pedigrees from the data
    EXCLUDE_PEDIGREES                  
    # Ignore 0 or more SNPs (1...N) from the data
    EXCLUDE_LOCUS                      
    # Maximum amount of missing data before a SNP is ignored from analysis
    MISSING_THRESHOLD                  0.1


    ########################## Basic Settings ###################
    # Describes the type of models of interest. Models consist of 1 or more SNPs.
    # The minimum number of SNPs in a model to be investigate.
    # Valid Settings: 1..COMBO_END
    COMBO_START                        1
    # The maxmimum number of SNPs to be considered.
    # Valid Settings: [COMBO_START..MAX_INT)
    COMBO_STOP                         2
    # Set the number of cross validation folds are used in analysis
    # Recommend settings: 1, 5, 10 (1 is no cross validation)
    CROSSVALINTERVAL                   5
    # Maximum number of models to be reported
    REPORTMODELCOUNT                   5
    # Threshold for reporting models. At most REPORTMODELCOUNT will be reported (those are sorted according to the T-Statistic
    REPORT_THRESHOLD                   0
    # Write out details regarding the contents of each cross validation fold
    VERBOSE_FOLDING                    0


    ######################### Permutation Tests ##################
    # Number of permutation runs to be executed. 1000 is recommended
    PTEST_COUNT                        1000
    # The seed associated with the tests. Each test gets a new seed
    PTEST_SEED                         1397
    # Short circuit the permutation tests. See manual for explaination
    PTEST_SHORTCIRCUIT                 0
    # How many simultaneous threads will be run
    # Each PTest can theoretically be run in it's own thread (Multiple threads
    # will not benefit if no ptests are being performed)
    THREAD_COUNT                       ${Runtime.runtime.availableProcessors()}
    # Determine how to handle mendelian errors when encountered. 
    # Acceptable Values:
    #   1 - Report errors, but do nothing
    #   2 - Report errors, and zero out loci in families where genotyping error has been found
    #   3 - Report errors and remove pedigrees where the number of genotyping errors exceeds threshold
    MENDELIAN_ERROR_LEVEL              1
    # Set the threshold, if level is 3
    MENDELIAN_PEDIGREE_THRESHOLD       0


    ######################### Report Names ########################
    # Pedigree report (genotypes and folding details)
    EXT_PEDIGREE                       pedigree
    # Distribution report (each item from the distribution)
    EXT_DISTRIBUTION                   dist
    # Write a copy of the 'cleaned' output to file
    WRITE_CLEAN_DATAFILE               No
    # Extension of the file (above)
    EXT_CLEANED_DATA                   ped
    EOF

    mdrpdt64 parameters.txt
    """

}
