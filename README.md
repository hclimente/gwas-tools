# gwas-tools

Short pipelines to deal with GWAS data.

# Impute

```
impute --ped data/example.ped --map dataexample.map --strand_info data/strand_info.txt --population EUR -with-docker hclimente/gwas-tools
```

# Run VEGAS 2ls

```
run_vegas --ped data/example.ped --map dataexample.map -with-docker hclimente/gwas-tools
```

# Map SNPs to genes from GENCODE

```
snp2gene --map data/example.map --genome GRCh38 -with-docker hclimente/gwas-tools
```

Dependencies:
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)
