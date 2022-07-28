#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")
bfile = tuple(bed, bim, fam)

params.covars = ''

process chisquare_test {

    input:
        tuple file(BED), file(BIM), file(FAM)

    output:
        file 'snp_association.tsv'
    
    """
    plink --bed ${BED} --bim ${BIM} --fam ${FAM} --assoc --allow-no-sex
    awk 'NR > 1 && \$9 != "NA" { print \$2,\$9 }' OFS='\\t' plink.assoc >snp_association.tsv
    """

}

process adjusted_logistic_regression {

    input:
        tuple file(BED), file(BIM), file(FAM)
        file COVARS

    output:
        file 'snp_association.tsv'
    
    """
    plink --bed ${BED} --bim ${BIM} --fam ${FAM} --logistic --covar ${COVAR}
    awk 'NR > 1 && \$5 == "ADD" && \$9 != "NA" { print \$2,\$9 }' OFS='\\t' plink.assoc.logistic >snp_association.tsv
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
