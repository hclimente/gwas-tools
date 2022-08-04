plink -bfile ${1} -recode oxford

cut -d' ' -f2,6- plink.gen | sed 's/ /\t/g' >X.tsv
cut -d' ' -f5 plink.sample >y.tsv

rm plink.gen plink.sample