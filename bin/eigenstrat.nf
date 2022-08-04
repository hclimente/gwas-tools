#!/usr/bin/env nextflow

params.out = '.'

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

process correct_population_structure {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file BED from bed
        file bim
        file fam

    output:
        file 'eigenstrat.plot.pdf'
        file 'eigenstrat.lambda'

    """
    plink --bfile ${BED.baseName} --recode --out ped
    smartpca.perl -i ped.ped -a ped.map -b ped.ped -o eigenstrat.pca -p eigenstrat.plot -e eigenstrat.eval -l eigenstrat.pca.log
    smarteigenstrat.perl -i ped.ped -a ped.map -b ped.ped -p eigenstrat.pca -l eigenstrat.eigenstrat.log -o eigenstrat.chisq
    gc.perl eigenstrat.chisq eigenstrat.lambda
    """

}
