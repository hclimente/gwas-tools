#!/usr/bin/env nextflow

params.out = "."

params.covar = ''
params.pheno = ''

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")

// parameters
params.N = 10

if (params.pheno != '') {

     pheno = file(params.pheno)
     original_fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

     process replace_phenotype {

         input:
             file PHENO from pheno
             file FAM from original_fam
             val I from params.i

         output:
             file 'fam' into fam

         script:
         template 'io/replace_phenotype.sh'

     }

} else {
    fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")
}

if (params.covar == '') {

   process run_regression {

       input:
           each I from 1..params.N
           val N from params.N
           file BED from bed
           file BIM from bim
           file FAM from fam

       output:
           file "part.epi.??.${I}" into parts

       """
       plink --bfile ${BED.baseName} --epistasis --parallel ${I} ${N} --allow-no-sex --out part
       """

    }

} else {

    covar = file(params.covar)

   process run_regression_with_covars {

       input:
           each I from 1..params.N
           val N from params.N
           file BED from bed
           file BIM from bim
           file FAM from fam
           file COVAR from covar

       output:
           file "part.epi.??.${I}" into parts

       """
       plink --bfile ${BED.baseName} --epistasis --covar ${COVAR} --parallel ${I} ${N} --allow-no-sex --out part
       """

    }


}

process join_parts {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file "*" from parts.collect()

    output:
        file "scored_interactions.regression.txt"

    """
    cat `ls -v part.epi.*` | sed 's/^ \\+//' | sed 's/ \\+/\t/g' | sed 's/\t\$//' >scored_interactions.regression.txt
    """

}
