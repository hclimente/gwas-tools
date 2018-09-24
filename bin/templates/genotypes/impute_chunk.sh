# SAMPLE, GEN, GENTIC_MAP, HAPS, LEGEND, STRAND_INFO, SAMPLES_REFERENCE
# imputed.gen, imputed.gen_warnings

chrXflags=''
if [[ ${CHR} == *X* ]]
  then
    chrXflags="-sample_g ${SAMPLE} -chrX"
    if [[ ${CHR} != *NONPAR* ]]
    then
        chrXflags="\$chrXflags -Xpar"
    fi
fi

grep "^${CHR} " ${STRAND_INFO} | cut -d' ' -f2,3 >strand_g_file

# create a gen with only the chromosome (the chromosomic region for efficiency)
# else, impute2 doesn't detect any snp
grep "^${CHR} " ${GEN} | awk '\$3 >= (${START} - 250000)'| awk '\$3 <= (${END} + 250000)' >region.gen
cut -f3 -d' ' region.gen | sort | uniq >pos_gen

zcat ${REFERENCE}/1000GP_Phase3_chr${CHR}.legend.gz | awk '\$2 >= (${START} - 250000)'| awk '\$2 <= (${END} + 250000)' >region.legend
cut -d' ' -f2 region.legend | sort | uniq >pos_leg

sort pos_* | uniq -d > repeated

if [ -s repeated ]
then
    # impute only on the screened SNPs
    ## only genotyped snps are included in the panel
    impute2 \
    -g region.gen \
    -m ${REFERENCE}/genetic_map_chr${CHR}_combined_b37.txt \
    -int ${START} ${END} \
    -h ${REFERENCE}/1000GP_Phase3_chr${CHR}.hap.gz \
    -l ${REFERENCE}/1000GP_Phase3_chr${CHR}.legend.gz \
    -Ne 20000 \
    -verbose \
    \$chrXflags \
    -k_hap ${SAMPLES_REFERENCE} \
    -os 2 -o imputed.gen \
    -strand_g strand_g_file
else
    mv region.gen imputed.gen
fi