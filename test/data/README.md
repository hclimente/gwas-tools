Toy example for testing purposes. It consists on a GWAS of a binary phenotype on 33 SNPs and 250 samples ([gwas.bed](gwas.bed), [gwas.bim](gwas.bim) and [gwas.fam](gwas.fam)), which map positionally to 18 genes ([snp2gene.tsv](snp2gene.tsv)). The 6 SNPs called "cX" (where X is an integer) are linked to the phenotype (see [assoc.chisq](assoc.chisq)); the remaining 27 ("oX") are not.

The causal SNPs map to the following genes: *ADM2*, *EP300*, *MAPK1*, *CHEK2*, *FBLN1* and *RBX1* (GENCODE 41, GRCh38). Hence, these genes also have high association scores with the phenotype([vegas2.tsv](vegas2.tsv)). Additionally, these genes form a clique in a gene-gene network ([edgelist.tsv](edgelist.tsv)).