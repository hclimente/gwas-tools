#!/usr/bin/env nextflow

params.out = '.'
params.pheno = ''

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")

// parameters
params.N = 1000
params.mode = 'binary' // binary or continuous

if (params.pheno != '') {

    pheno = file(params.pheno)
    original_fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

    process replace_phenotype {

        input:
            file PHENO from pheno
            file FAM from original_fam
            val I from params.i

        output:
            file "${FAM}" into fam

        script:
        template 'io/replace_phenotype.sh'

    }

} else {
    fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")
}

process bed2ped {

    input:
        file BED from bed
        file BIM from bim
        file FAM from fam

    output:
        file 'out.ped' into ped
        file 'out.map' into map

    script:
    template 'io/bed2ped.sh'

}

process ped2mbmdr {

    input:
        file PED from ped
        file MAP from map
        val MODE from params.mode

    output:
        file "mbmdrFile" into mbmdrFile
        file "labels" into mbmdrLabels

    """
    echo -e 'chr\tname\tcM\tpos' >fixed.map
    cat ${MAP} >>fixed.map

    mbmdr --plink2mbmdr --${MODE} -ped ${PED} -map fixed.map -o mbmdrFile -tr labels
    """

}

mbmdrFile.into { mbmdrFile_1; mbmdrFile_2; mbmdrFile_3; mbmdrFile_4 }

process compute_top_vectors {

    errorStrategy 'retry'
    maxRetries 3

    input:
        each I from 1..params.N
        val N from params.N
        file INPUT from mbmdrFile_1
        val MODE from params.mode

    output:
        file "top*.txt" into topFiles

    """
    mbmdr --${MODE} --gammastep1 -r 2019 -i ${I} -N ${N} ${INPUT}
    """

}

process merge_top_vectors {

    input:
        file "top*.txt" from topFiles.collect()
        file INPUT from mbmdrFile_2
        val MODE from params.mode
        val N from params.N

    output:
        file "topFile.txt" into topFile

    """
    mbmdr --${MODE} --gammastep2 -r 2019 -N ${N} ${INPUT}
    """

}

topFile.into { topFile_permutations; topFile_4 }

process compute_permutations {

    errorStrategy 'retry'
    maxRetries 3

    input:
        each I from 1..params.N
        file TOP from topFile_permutations
        file INPUT from mbmdrFile_3
        val MODE from params.mode

    output:
        file "perm_*" into permutFiles

    """
    mbmdr --${MODE} --gammastep3 -r 2019 -p 1 -o perm_${I}.txt -t ${TOP} ${INPUT}
    """

}

process merge_permutations {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file "perm_*.txt" from permutFiles.collect()
        val N from params.N
        file TOP from topFile_4
        file INPUT from mbmdrFile_4
        val MODE from params.mode

    output:
        file "scored_interactions.mbmdr.txt"

    """
    mbmdr --${MODE} --gammastep4 -r 2019 -c perm_ -q ${N} -t ${TOP} ${INPUT}
    mv *_output.txt scored_interactions.mbmdr.txt
    """

}
