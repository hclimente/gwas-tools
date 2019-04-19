# GEN, SAMPLE
# out.bed, out.bim, out.fam

plink -gen ${GEN} -sample ${SAMPLE} -hard-call-threshold 0.05 -make-bed -out out
