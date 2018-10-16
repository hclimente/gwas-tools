# gwas-tools

Containerized pipelines to deal with GWAS data.

# Installation

The easiest way to install gwas-tools is cloning the repository, and adding the bin folder to your path:

```
git clone git@github.com:hclimente/gwas-tools.git
export PATH=$PATH:$PWD/gwas-tools/bin
```

## Dependencies:

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)

# Impute

```
impute --ped data/example.ped --map dataexample.map --strand_info data/strand_info.txt --population EUR -with-docker hclimente/gwas-tools
```

# Run VEGAS 2

```
run_vegas --ped data/example.ped --map dataexample.map -with-docker hclimente/gwas-tools
```

# Map SNPs to genes from GENCODE

```
snp2gene --map data/example.map --genome GRCh38 -with-docker hclimente/gwas-tools
```


