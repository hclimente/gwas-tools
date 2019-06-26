# FAM, PHENO, I
# FAM

mv ${FAM} fam

cut -f1-5 -d' ' fam >ids
cut -f${I + 2} ${PHENO} >pheno

paste ids pheno >${FAM}
