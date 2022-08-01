# FAM, K, I
# samples_*.txt

get_seeded_random()
{
    seed="\$1"
    openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt \
    </dev/zero 2>/dev/null
}

total_samples=\$(grep -c '^' ${FAM})
n_samples=\$(expr \$total_samples - \$total_samples / $K)
cut -f1,2 -d' ' ${FAM} | shuf --random-source=<(get_seeded_random ${I}) | tail -n \$n_samples >samples_${I}.txt

plink -bfile ${BED.baseName} -keep samples_${I}.txt -allow-no-sex -make-bed -out train_${I}
plink -bfile ${BED.baseName} -remove samples_${I}.txt -allow-no-sex -make-bed -out test_${I}

