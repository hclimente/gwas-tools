# gwas-tools

*This repository is a fork from [hclimente/gwas-tools](https://github.com/hclimente/gwas-tools) with for now only the GWAS analysis part to prepare its adaptation for another project from Chloé-Agathe Azencott's team](https://cazencott.info/) at the [CBIO](https://cbio.ensmp.fr/)*

gwas-tools contains pipelines for common use-cases when dealing with GWAS datasets, from data preprocessing to biomarker discovery. 

## Installation

Start with cloning the repository, and adding the bin folder to your path:

```
git clone git@github.com:kumquatum/gwas-tools.git
export PATH=$PATH:$PWD/gwas-tools/bin
```

Then, install all tools described in [Dependencies](#dependencies) or build your own docker image based on the `Dockerfile` provided (some tools are under copyright and prevent us from providing a docker image). 

```
# Being in gwas-tools folder
docker build -t <name_of_your_image> .
```

The docker image can then used in nextflow by adding the parameter `-with-docker <name_of_your_image>`.


## Test files

A partial minimal set of files is available in `test/data` to demonstrate the use of gwas-tools. For the SConES tool to function, the PPI file need to be downloaded and prepared as in `bin/templates/dbs/biogrid.sh`

# Functions

- GWAS
    - [Data preprocessing](#data_preprocessing)
    - [Network-guided GWAS](#network_gwas)
<!--- Epistasis detection-->
<!--    - [Network-guided epistasis detection](#network_epistasis)-->


## GWAS

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

The available methods are:

- dmGWAS: 
```
dmgwas.nf --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 -with-docker <name_of_your_image>
```
- heinz: 
```heinz.nf --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 --fdr 0.5 -with-docker <name_of_your_image>
```
- HotNet2: 
```hotnet2.nf --scores test/data/scored_genes.vegas.txt --tab2 test/data/tab2 --hotnet2_path hotnet2 --lfdr_cutoff 0.125 -with-docker <name_of_your_image>
```
- LEAN: 
```lean.nf --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 -with-docker <name_of_your_image>
```
- SConES: 
```
old_scones.nf --bfile test/data/example --network gi --snp2gene test/data/snp2gene.tsv --tab2 test/data/tab2 -with-docker <name_of_your_image>
```
- Sigmod: 
```
sigmod.nf --sigmod SigMod_v2 --vegas test/data/scored_genes.vegas.txt --tab2 test/data/tab2 -with-docker <name_of_your_image>
```

<!--## Epistasis detection-->

<!--### Network-guided epistasis detection-->
<!--<a name="network_epistasis"></a>-->

<!--We developed a modular method to discover epistasis along the edges of a biological network. This facilitates the interpretability of the results and allows to find interactions that would otherwise be overcome by the multiple test burden. If you use this algorithm, please cite the following article:-->

<!--> Duroux D, Climente-González H, Azencott C-A, Van Steen K (2021) **Interpretable network-guided epistasis detection.** In press. https://doi.org/10.1101/2020.09.24.310136-->

<!--Usage:-->

<!--```-->
<!--network_epistasis.nf --bfile ../test/data/example --tab2 ../test/data/tab2 --snp2gene ../test/data/snp2gene.tsv --nperm 10-->
<!--```-->

# Dependencies

## To run the pipelines

* Nextflow
* Docker (optionnal)

## Pipeline itself

| Tool                                                                                                          | License                                                                   |
|---------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
<!--| [AntEpiSeeker](http://nce.ads.uga.edu/~romdhane/AntEpiSeeker/)                                                | ?                                                                         |-->
<!--| BEAM                                                                                                          | ?                                                                         |-->
| BEDOPS                                                                                                        | [GPLv2](https://github.com/bedops/bedops/blob/master/LICENSE)             |
<!--| [Biofilter](https://ritchielab.org/research/research-areas/expert-knowledge-bioinformatics/methods/biofilter) | ?                                                                         |-->
<!--| [GTOOL](https://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html)                                         | Copyright                                                                 |-->
| [HotNet2](https://github.com/raphael-group/hotnet2)                                                           | [Copyright](https://github.com/raphael-group/hotnet2/blob/master/LICENSE) |
| IMPUTE                                                                                                        | Copyright                                                                 |
<!--| liftOver                                                                                                      | [Copyright](http://hgdownload.soe.ucsc.edu/admin/exe/)                    |-->
<!--| [MB-MDR](http://bio3.giga.ulg.ac.be/index.php/software/mb-mdr/)                                               | ?                                                                         |-->
| [PLINK 1.90](https://www.cog-genomics.org/plink/1.9)                                                          | [GPLv3](https://www.cog-genomics.org/plink/1.9/general_usage)             |
<!--| [SMMB](https://www.ls2n.fr/listelogicielsequipe/DUKe/128/)                                                    | Freeware                                                                  |-->
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
