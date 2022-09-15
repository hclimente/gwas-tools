Toy example for testing purposes. It consists on a small GWAS (see gwas.bed, gwas.bim and gwas.fam) on 33 SNPs, which map positionally to 18 genes (see snp2gene.tsv). The 6 SNPs called "cX" (where X is an integer) linked to the phenotype (see assoc.chisq); the remaining 27 ("oX") are not.

The causal SNPs map to the following genes: *ADM2*, *EP300*, *MAPK1*, *CHEK2*, *FBLN1* and *RBX1* (GENCODE 41, GRCh38). Hence, these genes also have high association scores with the phenotype(see vegas2.tsv). Additionally, these genes form a clique in a gene-gene network (see edgelist.tsv).
