# gwas-tools

Pipelines for common use-cases when dealing with GWAS datasets, from preprocessing to biomarker discovery. All dependencies are provided in a [homonymous Docker image](https://hub.docker.com/repository/docker/hclimente/gwas-tools).

## Installation

Simply clone the repository and add the `bin/` folder to your path:

```bash
git clone git@github.com:hclimente/gwas-tools.git
export PATH=$PATH:$PWD/gwas-tools/bin
```

Install also the two dependencies:

- [Nextflow](https://www.nextflow.io/)
- Choose one of:
    - [Docker](https://docs.docker.com/desktop/#download-and-install)
    - [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

Each pipeline uses multiple tools and libraries, which I compile in a Docker image. Those that can be freely distributed are available in [hclimente/gwas-tools](https://hub.docker.com/r/hclimente/gwas-tools). You will need to include the remaining ones into a new Docker image (hclimente/gwas-tools-extra) by running:

```bash
make docker
```

Then, simply add '*-with-docker hclimente/gwas-tools[-extra]*' or '*-with-singularity hclimente/gwas-tools[-extra]*' when launching a pipeline.

# Functions

- GWAS
    - [Data preprocessing](#data_preprocessing)
    - [Network-guided GWAS](#network_gwas)
- Epistasis detection
    - [Network-guided epistasis detection](#network_epistasis)

## GWAS

### Data preprocessing
<a name="data_preprocessing"></a>

- Impute a dataset: `impute --bfile test/data/example --strand_info test/data/strand_info.txt --population EUR -with-docker hclimente/gwas-tools`
- Run VEGAS2: `vegas2.nf --bfile test/data/example --gencode 31 --genome 37 --buffer 50000 --vegas_params '-top 10' -with-docker hclimente/gwas-tools`
- Map SNPs to GENCODE genes: `snp2gene.nf --bim test/data/example.map --genome GRCh38 -with-docker hclimente/gwas-tools`

### Network-guided GWAS
<a name="network_gwas"></a>

We adapted and benchmarked multiple algorithms for the detection of SNPs associated to a phenotype. If you use any of the following algorithms, please cite the original algorithm's article, and the following article:

> Climente-González H, Lonjou C, Lesueur F, GENESIS study group, Stoppa-Lyonnet D, et al. (2021) **Boosting GWAS using biological networks: A study on susceptibility to familial breast cancer.** PLOS Computational Biology 17(3): e1008819. https://doi.org/10.1371/journal.pcbi.1008819

The available methods are:

- dmGWAS: `dmgwas.nf --vegas scored_genes.top10.txt --edgelist test/data/interactions.tsv -with-docker hclimente/gwas-tools`
- heinz: `heinz.nf --vegas scored_genes.top10.txt --edgelist test/data/interactions.tsv --fdr 0.5 -with-docker hclimente/gwas-tools`
- HotNet2: `hotnet2.nf --scores scored_genes.top10.txt --edgelist test/data/interactions.tsv --hotnet2_path hotnet2 --lfdr_cutoff 0.125 -with-docker hclimente/gwas-tools`
- LEAN: `lean.nf --vegas scored_genes.top10.txt --edgelist test/data/interactions.tsv -with-docker hclimente/gwas-tools`
- SConES: `scones.nf --bfile test/data/example --network gi --snp2gene test/data/snp2gene.tsv --edgelist test/data/interactions.tsv -with-docker hclimente/gwas-tools`
- Sigmod: `sigmod.nf --sigmod SigMod_v2 --vegas scored_genes.top10.txt --edgelist test/data/interactions.tsv -with-docker hclimente/gwas-tools`

## Epistasis detection

### Network-guided epistasis detection
<a name="network_epistasis"></a>

We developed a modular method to discover epistasis along the edges of a biological network. This facilitates the interpretability of the results and allows to find interactions that would otherwise be overcome by the multiple test burden. If you use this algorithm, please cite the following article:

> Duroux D, Climente-González H, Azencott C-A, Van Steen K (2021) **Interpretable network-guided epistasis detection.** In press. https://doi.org/10.1101/2020.09.24.310136

Usage:

```bash
network_epistasis.nf --bfile ../test/data/example --tab2 ../test/data/interactions.tab2 --snp2gene ../test/data/snp2gene.tsv --nperm 10
```
