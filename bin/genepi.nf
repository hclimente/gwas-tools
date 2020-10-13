#!/usr/bin/env nextflow

params.out = "."
pheno_type = (params.phenotype == 'discrete')? 'c' : 'r'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

process bed2gen {

    input:
    file BED from bed
    file BIM from bim
    file FAM from fam

    output:
    file 'out.gen' into gen
    file 'out.sample' into sample

    script:
    template 'io/bed2gen.sh'

}

process make_pheno {

    input:
    file FAM from fam

    output:
    file 'pheno.csv' into pheno

    script:
    """
    cut -d' ' -f6 ${FAM} >pheno.csv
    """

}

process genepi {

    publishDir "$params.out", overwrite: true, mode: 'copy'

    input:
    file GEN from gen
    file PHENO from pheno
    val PHENO_TYPE from pheno_type 

    output:
    file 'crossGeneResult/Result.csv'

    script:
    """
    GenEpi -g ${GEN} -p ${PHENO} -m ${PHENO_TYPE} -o ./
    """

}
