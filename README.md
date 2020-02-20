# gwas-tools

[![Docker Pulls](https://img.shields.io/docker/cloud/build/hclimente/gwas-tools.svg?style=popout-square&logo=docker)](https://cloud.docker.com/swarm/hclimente/repository/docker/hclimente/gwas-tools)

This repository contains pipelines for common use-cases when dealing with GWAS datasets, from data preprocessing to biomarker discovery. 

## Installation

The easiest way to install gwas-tools is cloning the repository, and adding the bin folder to your path:

```
git clone git@github.com:hclimente/gwas-tools.git
export PATH=$PATH:$PWD/gwas-tools/bin
```

The pipelines are written in [Nextflow](https://www.nextflow.io/), and makes use of multiple tools (see [Dependencies](#dependencies)). These tools need to be installed independently on a per-pipeline basis. However, those that can be distributed are included in a Docker image, which can be used adding the parameter '*-with-docker hclimente/gwas-tools*' or '*-with-singularity hclimente/gwas-tools*'.

# Functions

- Data preprocessing
    - [Impute a dataset](#impute-a-dataset)
    - [PCA on a dataset]()
    - [Run VEGAS2 on a dataset](#run-vegas2)
    - [Map SNPs to genes](#map-snps-to-gencode-genes)
    - [Map eQTLs to their regulated genes]()
    - [Filter SNPs with Biofilter]()
    - [Lift coordinates]()

- Network GWAS
    - [dmGWAS]()
    - [heinz]()
    - [HotNet2]()
    - [LEAN]()
    - [SConES]()
    - [Sigmod]()

- Epistasis detection
    - [MB-MDR]()

## Data preprocessing

### Impute a dataset

```
impute --bfile test/data/example --strand_info test/data/strand_info.txt --population EUR -with-docker hclimente/gwas-tools
```

### Run VEGAS2

```
vegas2.nf --bfile test/data/example -with-docker hclimente/gwas-tools
```

### Map SNPs to GENCODE genes

```
snp2gene.nf --bim test/data/example.map --genome GRCh38 -with-docker hclimente/gwas-tools
```

## Network GWAS

# Dependencies

| Tool         | Docker'd | License   
| -------------|:---------|:----------
| AntEpiSeeker | No       | ?         
| BEAM         | No       | ?         
| BEDOPS       | Yes      | [GPLv2](https://github.com/bedops/bedops/blob/master/LICENSE)
| [Biofilter](https://ritchielab.org/research/research-areas/expert-knowledge-bioinformatics/methods/biofilter)    | No       | ?
| [GTOOL](https://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html) | No       | Copyright
| [HotNet2](https://github.com/raphael-group/hotnet2) | No       | [Copyright](https://github.com/raphael-group/hotnet2/blob/master/LICENSE)
| IMPUTE       | No       | Copyright
| liftOver     | No       | [Copyright](http://hgdownload.soe.ucsc.edu/admin/exe/)
| [MB-MDR](http://bio3.giga.ulg.ac.be/index.php/software/mb-mdr/) | No       | ?
| [PLINK 1.90](https://www.cog-genomics.org/plink/1.9) | Yes      | [GPLv3](https://www.cog-genomics.org/plink/1.9/general_usage)
| R::BioNet    | Yes      | [GPLv2](https://bioconductor.org/packages/release/bioc/html/BioNet.html)
| R::dmGWASv3  | Yes      | [GPLv2](https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0-manual.pdf)
| R::LEANR     | Yes      | [GPLv3](https://cran.r-project.org/web/packages/LEANR/)
| R::martini   | Yes      | [MIT](https://bioconductor.org/packages/release/bioc/html/martini.html)
| R::ranger    | Yes      | [GPLv3](https://cran.r-project.org/web/packages/ranger/)
| R::SigModv2  | No       | ?
| [VEGAS2v02](https://vegas2.qimrberghofer.edu.au/) | No       | ?
