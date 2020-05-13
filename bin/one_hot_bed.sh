plink -bfile ${1} -recode oxford

cut -d' ' -f2,6- plink.gen >X.dat
cut -d' ' -f5 plink.sample >y.dat