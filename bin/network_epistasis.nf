#!/usr/bin/env nextflow

params.out = "."
params.nperm = 1000

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

tab2 = file("${params.tab2}")
snp2gene = file("${params.snp2gene}")

process preprocess_data {

	input:
        file BED from bed
        file bim
        file FAM from fam
	file TAB2 from tab2 
	file SNP2GENE from snp2gene
	val I from params.nperm

	output:
	set 'post_qc.bed', 'post_qc.bim', 'post_qc.fam', 'pheno' into filtered_bed
	file 'studied_snp_pairs.tsv' into studied_snp_pairs

	"""
	# TODO DOUBLE CHECK create SNP pairs
	R -e 'library(tidyverse); snp2gene <- read_tsv("${SNP2GENE}", col_types = "cc"); read_tsv("${TAB2}", col_types = cols(.default = col_character())) %>% inner_join(snp2gene, by = c("Official Symbol Interactor A" = "gene")) %>% inner_join(snp2gene, by = c("Official Symbol Interactor B" = "gene"), suffix = c("_1", "_2")) %>% select(snp_1, snp_2) %>% write_tsv("studied_snp_pairs.tsv")'

	# QC + keep only SNPs in models
	cut -f1 studied_snp_pairs.tsv >tmp
	cut -f2 studied_snp_pairs.tsv >>tmp
	sort tmp | uniq >included_snps
	plink -bfile ${BED.baseName} -extract included_snps -maf 0.05 -hwe 0.001 -make-bed -out filtered 

	# LD pruning
	plink -bfile filtered -indep-pairwise 50 5 0.75 -out snps_75
	plink -bfile filtered -extract snps_75.prune.in -make-bed -out post_qc

	# create pheno file
	cut -d' ' -f6 ${FAM} >orig
	cut -d' ' -f1,2,6 ${FAM} | sed 's/ /\t/g'  >pheno
	for i in {1..${I}}
	do
		shuf orig >tmp1
		paste pheno tmp1 >tmp2
		mv tmp2 pheno
	done
	"""

}

process snp_epistasis {
	
	input:
        each I from 1..(params.nperm + 1)
	set file(BED), file(BIM), file(FAM), file(PHENO) from filtered_bed 

	output:
        set val(I), 'scored_interactions.regression.txt' into snp_pairs, snp_pairs_null 

        """
        epistasis_regression.nf --bfile ${BED.baseName} --pheno ${PHENO} --i ${I} 
	""" 

}

snp_pairs_null
	.filter { it -> it[0] != 1 }
	.flatMap { it -> it[1] }
	.set { snp_pairs_null_filtered }

process gene_epistasis {

	input:
	file '*' from snp_pairs_null_filtered .collect()
	set I, file(SNP_PAIRS) from snp_pairs
	file SNP2GENE from snp2gene
	file SNP2SNP from studied_snp_pairs
        file TAB2 from tab2

	output:
	set val(I), 'scored_gene_pairs.tsv' into scored_gene_pairs, scored_gene_pairs_null

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)

	snp2snp <- read_tsv('$SNP2SNP', col_types = 'cc') %>%
        	mutate(uniq_snp_id = cbind(SNP_1, SNP_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
		select(uniq_snp_id)
	threshold <- lapply(list.files(pattern = 'permuted_association'), function(x) {
				read_tsv(x) %>%
					mutate(uniq_snp_id = cbind(SNP_1, SNP_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
					inner_join(snp2snp, by = 'uniq_snp_id') %>%
					arrange(P) %>%
					top_n(1) %>%
					select(P)
			   }) %>%
		do.call(bind_rows, .) %>%
		.$P %>%
		quantile(.95)

	snp_pairs <- read_tsv('${SNP_PAIRS}') %>%
        	mutate(uniq_snp_id = cbind(SNP_1, SNP_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
		filter(P < threshold) %>%
		inner_join(snp2snp, by = 'uniq_snp_id') %>%
		select(uniq_snp_id, P)

	read_tsv('${TAB2}', col_types = cols(.default = col_character())) %>%
		inner_join(snp2gene, by = c('Official Symbol Interactor A' = 'gene')) %>%
                inner_join(snp2gene, by = c('Official Symbol Interactor B' = 'gene')) %>%
                mutate(uniq_snp_id = cbind(SNP_1, SNP_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_'),
		       uniq_gene_id = cbind(gene_1, gene_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
		inner_join(snp_pairs, by = 'uniq_snp_id') %>%
		group_by(uniq_gene_id) %>%
		summarize(tau_0.05 = prod(P[P < 0.05]),
			  tau_0.01 = prod(P[P < 0.01]),
			  tau_0.001 = prod(P[P < 0.001])) %>%
		write_tsv('scored_gene_pairs.tsv')
	"""

}

// TODO create two subsets
scored_gene_pairs
        .filter { it -> it[0] == 1 }
        .flatMap { it -> it[1] }
        .set { gene_pairs }
scored_gene_pairs_null
        .filter { it -> it[0] != 1 }
        .flatMap { it -> it[1] }
        .set { gene_pairs_null }

process involved_pathways {

	input:
	file 'null_*' from gene_pairs_null .collect()
	file GENE_PAIRS from gene_pairs
	file TAB2 from tab2

	output:
	file 'scored_gene_pairs.tsv'
	file 'scored_pathways.tsv'

	"""
        #!/usr/bin/env Rscript

	library(tidyverse)
	library(igraph)

        net <- read_tsv('${TAB2}', col_types = cols(.default = col_character())) %>%
                select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
		graph_from_data_frame
	"""

}
