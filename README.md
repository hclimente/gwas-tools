# gwas-tools

[![Docker Pulls](https://img.shields.io/docker/cloud/build/hclimente/gwas-tools.svg?style=popout-square&logo=docker)](https://hub.docker.com/repository/docker/hclimente/gwas-tools)

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
    - [dmGWAS](#dmgwas)
    - [heinz](#heinz)
    - [HotNet2](#hotnet2)
    - [LEAN](#lean)
    - [SConES](#scones)
    - [Sigmod](#sigmod)

- Epistasis detection
    - [MB-MDR]()

## Data preprocessing

### Impute a dataset

```
impute --bfile test/data/example --strand_info test/data/strand_info.txt --population EUR -with-docker hclimente/gwas-tools
```

### Run VEGAS2

```
vegas2.nf --bfile test/data/example --gencode 31 --genome 37 --buffer 50000 --vegas_params '-top 10' -with-docker hclimente/gwas-tools
```

### Map SNPs to GENCODE genes

```
snp2gene.nf --bim test/data/example.map --genome GRCh38 -with-docker hclimente/gwas-tools
```

## Network GWAS

### dmGWAS

```
dmgwas.nf --vegas scored_genes.top10.txt --tab2 test/data/interactions.tab2 -with-docker hclimente/gwas-tools
```

### heinz

```
heinz.nf --vegas scored_genes.top10.txt --tab2 test/data/interactions.tab2 --fdr 0.5 -with-docker hclimente/gwas-tools
```

### HotNet2

```
hotnet2.nf --scores scored_genes.top10.txt --tab2 test/data/interactions.tab2 --hotnet2_path hotnet2 --lfdr_cutoff 0.125 -with-docker hclimente/gwas-tools
```

### LEAN

```
lean.nf --vegas scored_genes.top10.txt --tab2 test/data/interactions.tab2 -with-docker hclimente/gwas-tools
```

### SConES

```
old_scones.nf --bfile test/data/example.map --network gi --snp2gene test/data/snp2gene.tsv --tab2 test/data/interactions.tab2 -with-docker hclimente/gwas-tools
```

### Sigmod

```
sigmod.nf --sigmod SigMod_v2 --vegas scored_genes.top10.txt --tab2 test/data/interactions.tab2 -with-docker hclimente/gwas-tools
```

# Dependencies

| Tool         | Docker'd | License   
| -------------|:---------|:----------
| [AntEpiSeeker](http://nce.ads.uga.edu/~romdhane/AntEpiSeeker/) | No       | ?         
| BEAM         | No       | ?         
| BEDOPS       | Yes      | [GPLv2](https://github.com/bedops/bedops/blob/master/LICENSE)
| [Biofilter](https://ritchielab.org/research/research-areas/expert-knowledge-bioinformatics/methods/biofilter)    | No       | ?
| [GTOOL](https://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html) | No       | Copyright
| [HotNet2](https://github.com/raphael-group/hotnet2) | No       | [Copyright](https://github.com/raphael-group/hotnet2/blob/master/LICENSE)
| IMPUTE       | No       | Copyright
| liftOver     | No       | [Copyright](http://hgdownload.soe.ucsc.edu/admin/exe/)
| [MB-MDR](http://bio3.giga.ulg.ac.be/index.php/software/mb-mdr/) | No       | ?
| [PLINK 1.90](https://www.cog-genomics.org/plink/1.9) | Yes      | [GPLv3](https://www.cog-genomics.org/plink/1.9/general_usage)
| [SMMB](https://www.ls2n.fr/listelogicielsequipe/DUKe/128/) | No       | Freeware
| [VEGAS2v02](https://vegas2.qimrberghofer.edu.au/) | Yes       | [GPLv3](https://vegas2.qimrberghofer.edu.au/vegas2v2)
| R::biglasso  | Yes      | [GPLv3](https://cran.r-project.org/web/packages/biglasso/)
| R::bigmemory | Yes      | [LGPLv3](https://cran.r-project.org/web/packages/bigmemory/)
| R::CASMAP    | Yes      | [GPLv2](https://cran.r-project.org/web/packages/CASMAP/index.html)
| R::BioNet    | Yes      | [GPLv2](https://bioconductor.org/packages/release/bioc/html/BioNet.html)
| R::dmGWASv3  | Yes      | [GPLv2](https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0-manual.pdf)
| R::igraph    | Yes      | [GPLv3](https://cran.r-project.org/web/packages/igraph/)
| R::LEANR     | Yes      | [GPLv3](https://cran.r-project.org/web/packages/LEANR/)
| R::martini   | Yes      | [MIT](https://bioconductor.org/packages/release/bioc/html/martini.html)
| R::ranger    | Yes      | [GPLv3](https://cran.r-project.org/web/packages/ranger/)
| R::SigModv2  | No       | ?
| R::SKAT      | Yes      | [GPLv3](https://cran.r-project.org/web/packages/SKAT/)
| R::snpStats  | Yes      | [GPLv3](http://bioconductor.org/packages/release/bioc/html/snpStats.html)
