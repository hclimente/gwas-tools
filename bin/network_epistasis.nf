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
	# TODO REMOVE --N 3
	epistasis_regression.nf --bfile ${BED.baseName} --pheno ${PHENO} --i ${I} --N 3 
	""" 

}

snp_pairs_null
	.filter { it -> it[0] != 1 }
	.flatMap { it -> it[1] }
	.set { snp_pairs_null_filtered }

process gene_epistasis {

	input:
		file 'permuted_association_*' from snp_pairs_null_filtered .collect()
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
		mutate(uniq_snp_id = cbind(snp_1, snp_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
		select(uniq_snp_id)
		snp2gene <- read_tsv('$SNP2GENE', col_types = 'cc')
	threshold <- lapply(list.files(pattern = 'permuted_association_'), function(x) {
			read_tsv(x, col_types = 'ccccddd') %>%
				mutate(uniq_snp_id = cbind(SNP1, SNP2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
				inner_join(snp2snp, by = 'uniq_snp_id') %>%
				arrange(P) %>%
				top_n(1) %>%
				select(P)
			}) %>%
		do.call(bind_rows, .) %>%
		.\$P %>%
		quantile(.95)

	snp_pairs <- read_tsv('${SNP_PAIRS}', col_types = 'ccccddd') %>%
		mutate(uniq_snp_id = cbind(SNP1, SNP2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
		filter(P < threshold) %>%
		inner_join(snp2snp, by = 'uniq_snp_id') %>%
		select(uniq_snp_id, P)

	read_tsv('${TAB2}', col_types = cols(.default = col_character())) %>%
		inner_join(snp2gene, ., by = c('gene' = 'Official Symbol Interactor A')) %>%
		inner_join(snp2gene, ., by = c('gene' = 'Official Symbol Interactor B'), suffix = c('_1','_2')) %>%
		rename(gene_1 = gene) %>%
		mutate(uniq_snp_id = cbind(snp_1, snp_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_'),
			   uniq_gene_id = cbind(gene_1, gene_2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
		inner_join(snp_pairs, by = 'uniq_snp_id') %>%
		group_by(uniq_gene_id) %>%
		summarize(tau_05 = prod(P[P < 0.05]),
				  tau_01 = prod(P[P < 0.01]),
				  tau_001 = prod(P[P < 0.001])) %>%
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

process pathway_epistasis {

	publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file 'permuted_association_*' from gene_pairs_null .collect()
		file GENE_PAIRS from gene_pairs
		file TAB2 from tab2

	output:
		file 'scored_gene_pairs.tsv'
		file 'scored_pathways.tsv'

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)
	library(igraph)

	# compute gene-pair association
	gene_pairs <- read_tsv('${GENE_PAIRS}') %>%
	mutate(experiment = 'original')
	studied_gene_pairs <- gene_pairs\$uniq_gene_id		

	gene_pairs_assoc <- lapply(list.files(pattern = 'permuted_association_'), function(x) {
			read_tsv(x) %>%
				filter(uniq_gene_id %in% studied_gene_pairs)
		}) %>%
		do.call(bind_rows, .) %>%
		mutate(experiment = 'permutation') %>%
		bind_rows(gene_pairs) %>%
		group_by(uniq_gene_id) %>%
		mutate(P_tau_05 = rank(tau_05) / n(),
			   P_tau_01 = rank(tau_01) / n(),
			   P_tau_001 = rank(tau_001) / n(),
			   min_P = pmin(P_tau_05, P_tau_01, P_tau_001),
			   P = rank(min_P) / n() ) %>%
		filter(experiment == 'original') %>%
		select(-experiment)
	write_tsv(gene_pairs_assoc, 'scored_gene_pairs.tsv')
		
	# compute pathway association
	net <- read_tsv('${TAB2}', col_types = cols(.default = col_character())) %>%
		select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
		graph_from_data_frame
 
	# TODO create background
	bg <- NULL

	lapply(filter(gene_pairs_assoc, P < 0.05)$uniq_gene_id, function(uniq_gene_id) {
			genes <- unlist(strsplit('a_b', split = '_'))
			# TODO more than 1?
			delete_edges(gsub('_', '|', 'uniq_gene_id', fixed = TRUE)) %>%
				shortest_paths(from = genes[1], to = genes[2]) %>%
				goenrich(bg = bg)
	}) %>%
	do.call(bind_rows, .) %>%
	write_tsv('scored_pathways.tsv')	
	"""

}