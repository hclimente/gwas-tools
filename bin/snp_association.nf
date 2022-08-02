#!/usr/bin/env nextflow

include { get_bfile } from './templates/utils.nf'

// gwas
bfile = get_bfile(params.bfile)

params.covars = ''

process chisquare_test {
    
    tag { BED.getBaseName() }

    input:
        tuple file(BED), file(BIM), file(FAM)

    output:
        path "${BED.getBaseName()}.chisq"
    
    """
    plink --bed ${BED} --bim ${BIM} --fam ${FAM} --assoc --allow-no-sex
    awk 'NR > 1 && \$9 != "NA" { print \$2,\$9 }' OFS='\\t' plink.assoc >${BED.getBaseName()}.chisq
    """

}

process adjusted_logistic_regression {

    tag { BED.getBaseName() }

    input:
        tuple path(BED), path(BIM), path(FAM)
        file COVARS

    output:
        path "${BED.getBaseName()}.logreg"
    
    """
    plink --bed ${BED} --bim ${BIM} --fam ${FAM} --logistic --covar ${COVAR}
    awk 'NR > 1 && \$5 == "ADD" && \$9 != "NA" { print \$2,\$9 }' OFS='\\t' plink.assoc.logistic >${BED.getBaseName()}.logreg
    """

}


workflow snp_association_nf {
    take:
        bfile
        covars
    main:
        if (covars == '') {
            chisquare_test(bfile)
            assoc = chisquare_test.out
        } else {
            adjusted_logistic_regression(bfile, file(covars))
            assoc = adjusted_logistic_regression.out
        }

    emit:
        assoc
}

workflow {
    main:
        snp_association_nf(bfile, params.covars)
    emit:
        snp_association_nf.out
}
