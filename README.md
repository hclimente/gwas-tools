# gwas-tools

[![Docker Pulls](https://img.shields.io/docker/cloud/build/hclimente/gwas-tools.svg?style=popout-square&logo=docker)](https://cloud.docker.com/swarm/hclimente/repository/docker/hclimente/gwas-tools)

This repository contains pipelines for common use-cases when dealing with GWAS datasets, from data preprocessing to biomarker discovery. 

## Installation

The easiest way to install gwas-tools is cloning the repository, and adding the bin folder to your path:

```
git clone git@github.com:hclimente/gwas-tools.git
export PATH=$PATH:$PWD/gwas-tools/bin
```

The pipelines are written in [Nextflow](https://www.nextflow.io/), and makes use of several tools. These tools need to be installed independently on a per-pipeline basis, or used from a Docker image using the option '*-with-docker hclimente/gwas-tools*' or '*-with-singularity hclimente/gwas-tools*'.

# Functions

- Data preprocessing
    - [Impute a dataset](#impute-a-dataset)
    - [PCA on a dataset]()
    - [Run VEGAS2 on a dataset](#run-vegas2)
    - [Map SNPs to genes](#map-snps-to-gencode-genes)
    - [Map eQTLs to their regulated genes]()
    - [Filter SNPs with Biofilter]()
    - [Lift coordinates]()

- Network-based discovery
    - [SConES]()
    - [dmGWAS]()
    - [LEAN]()
    - [Sigmod]()

- Epistasis detection
    - [MB-MDR]()

# Data preprocessing

## Impute a dataset

```
impute --bfile test/data/example --strand_info test/data/strand_info.txt --population EUR -with-docker hclimente/gwas-tools
```

## Run VEGAS2

```
run_vegas --bfile test/data/example -with-docker hclimente/gwas-tools
```

## Map SNPs to GENCODE genes

```
snp2gene --bim test/data/example.map --genome GRCh38 -with-docker hclimente/gwas-tools
```

# Dependencies

| Tool         | Docker'd | License   | Source                                                                          
| -------------|:---------|:----------|:--------------------
| AntEpiSeeker | No       | ?         |
| BEAM         | No       | ?         |
| BEDOPS       | Yes      | GPLv2     | [LICENSE](https://github.com/bedops/bedops/blob/master/LICENSE)
| Biofilter    | No       | ?
| GTOOL        | No       | Copyright | Downloaded license file
| HotNet2      | No       | Copyright | [LICENSE](https://github.com/raphael-group/hotnet2/blob/master/LICENSE)
| IMPUTE       | No       | Copyright | Downloaded license file
| liftOver     | No       | Copyright | [link](http://hgdownload.soe.ucsc.edu/admin/exe/)
| MB-MDR       | No       | ?
| PLINK 1.90   | Yes      | GPLv3     | [link](https://www.cog-genomics.org/plink/1.9/general_usage)
| R::BioNet    | Yes      | GPLv2     | [Bioconductor](https://bioconductor.org/packages/release/bioc/html/BioNet.html)
| R::dmGWASv3  | Yes      | GPLv2     | [link](https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0-manual.pdf)
| R::LEANR     | Yes      | GPLv3     | [CRAN](https://cran.r-project.org/web/packages/LEANR/)
| R::martini   | Yes      | MIT       | [Bioconductor](https://bioconductor.org/packages/release/bioc/html/martini.html)
| R::ranger    | Yes      | GPLv3     | [CRAN](https://cran.r-project.org/web/packages/ranger/)
| R::SigModv2  | No       | ?
| VEGAS2v02    | No       | ?
