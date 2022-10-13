def get_bfile(bfile) {

    bed = file("${bfile}.bed")
    bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
    fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")
    tuple(bed, bim, fam)

}

