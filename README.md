# gwas-tools

*This repository is a fork from [hclimente/gwas-tools](https://github.com/hclimente/gwas-tools) with for now only the GWAS analysis part to prepare its adaptation for another project from [Chloé-Agathe Azencott's team](https://cazencott.info/) at the [CBIO](https://cbio.ensmp.fr/)*

gwas-tools contains pipelines for common use-cases when dealing with GWAS datasets, from data preprocessing to biomarker discovery. 

  

- [Installation](#installation)
  * [gwas-tools core](#gwas-tools-core)
  * [Dependencies](#dependencies)
    + [To run the pipelines](#to-run-the-pipelines)
    + [Pipeline itself](#pipeline-itself)
  * [Test files](#test-files)
- [Functions](#functions)
  * [Data preprocessing](#data-preprocessing)
  * [Network-guided GWAS](#network-guided-gwas)
    + [Gene-based methods](#gene-based-methods)
    + [SNP based methods](#snp-based-methods)


## Installation

### gwas-tools core

Clone the repository, and add the bin folder to your path:
```
git clone git@github.com:kumquatum/gwas-tools.git
export PATH=$PATH:$PWD/gwas-tools/bin
```

### Dependencies

#### To run the pipelines

* Nextflow
* Docker (optionnal)

#### Pipeline itself

Install all tools described below or build your own docker image based on the `Dockerfile` provided (some tools are under copyright and prevent us from providing a docker image) with : 
```
# Being in gwas-tools folder
docker build -t <name_of_your_image> .
```
The docker image can then used in nextflow by adding the parameter `-with-docker <name_of_your_image>`.

| Tool                                                                                                          | License                                                                   |
|---------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| BEDOPS                                                                                                        | [GPLv2](https://github.com/bedops/bedops/blob/master/LICENSE)             |
| [HotNet2](https://github.com/raphael-group/hotnet2)                                                           | [Copyright](https://github.com/raphael-group/hotnet2/blob/master/LICENSE) |
| IMPUTE                                                                                                        | Copyright                                                                 |
| [PLINK 1.90](https://www.cog-genomics.org/plink/1.9)                                                          | [GPLv3](https://www.cog-genomics.org/plink/1.9/general_usage)             |
| [VEGAS2v02](https://vegas2.qimrberghofer.edu.au/)                                                             | [GPLv3](https://vegas2.qimrberghofer.edu.au/vegas2v2)                     |
| R::biglasso                                                                                                   | [GPLv3](https://cran.r-project.org/web/packages/biglasso/)                |
| R::bigmemory                                                                                                  | [LGPLv3](https://cran.r-project.org/web/packages/bigmemory/)              |
| R::CASMAP                                                                                                     | [GPLv2](https://cran.r-project.org/web/packages/CASMAP/index.html)        |
| R::BioNet                                                                                                     | [GPLv2](https://bioconductor.org/packages/release/bioc/html/BioNet.html)  |
| R::dmGWASv3                                                                                                   | [GPLv2](https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0-manual.pdf)             |
| R::igraph                                                                                                     | [GPLv3](https://cran.r-project.org/web/packages/igraph/)                  |
| R::LEANR                                                                                                      | [GPLv3](https://cran.r-project.org/web/packages/LEANR/)                   |
| R::martini                                                                                                    | [MIT](https://bioconductor.org/packages/release/bioc/html/martini.html)   |
| R::ranger                                                                                                     | [GPLv3](https://cran.r-project.org/web/packages/ranger/)                  |
| R::SigModv2                                                                                                   | ?                                                                         |
| R::SKAT                                                                                                       | [GPLv3](https://cran.r-project.org/web/packages/SKAT/)                    |
| R::snpStats                                                                                                   | [GPLv3](http://bioconductor.org/packages/release/bioc/html/snpStats.html) |

### Test files

A partial minimal set of files is available in `test/data` to demonstrate the use of gwas-tools. For the SConES tool to function, the PPI file need to be downloaded and prepared as in `bin/templates/dbs/biogrid.sh`


## Functions

### Data preprocessing
<a name="data_preprocessing"></a>

<!--- Impute a dataset: `impute --bfile test/data/example --strand_info test/data/strand_info.txt --population EUR -with-docker <name_of_your_image>`-->
- Run VEGAS2: 
```
vegas2.nf --bfile test/data/example --gencode 31 --genome 37 --buffer 50000 --vegas_params '-top 10' -with-docker <name_of_your_image>
```
- Map SNPs to GENCODE genes: 
```
snp2gene.nf --bim test/data/example.map --genome GRCh38 -with-docker <name_of_your_image>
```

### Network-guided GWAS
<a name="network_gwas"></a>

Multiple algorithms were adapted and benchmarked for the detection of SNPs associated to a phenotype. If you use any of the following algorithms, please cite the following article:

> Climente-González H, Lonjou C, Lesueur F, GENESIS study group, Stoppa-Lyonnet D, et al. (2021) **Boosting GWAS using biological networks: A study on susceptibility to familial breast cancer.** PLOS Computational Biology 17(3): e1008819. https://doi.org/10.1371/journal.pcbi.1008819

#### Gene-based methods 

- dmGWAS: 
```
dmgwas.nf --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 -with-docker <name_of_your_image>
```
- heinz: 
```
heinz.nf --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 --fdr 0.5 -with-docker <name_of_your_image>
```
- HotNet2: 
```
hotnet2.nf --scores test/data/scored_genes.vegas.txt --tab2 test/data/tab2 --hotnet2_path hotnet2 --lfdr_cutoff 0.125 -with-docker <name_of_your_image>
```
- LEAN: 
```
lean.nf --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 -with-docker <name_of_your_image>
```
- Sigmod: 
```
sigmod.nf --sigmod SigMod_v2 --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 -with-docker <name_of_your_image>
```

#### SNP based methods

- SConES: 
```
old_scones.nf --bfile test/data/example --network gi --snp2gene test/data/snp2gene.tsv --tab2 test/data/tab2 -with-docker <name_of_your_image>
```



