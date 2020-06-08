#!/usr/bin/env nextflow

params.out = "."
params.nperm = 999

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

// gene & snp networks
tab2 = file("${params.tab2}")
snp2gene = file("${params.snp2gene}")
excluded_snps = file("${params.excluded}")

process make_snp_models {

	input:
		file TAB2 from tab2 
		file SNP2GENE from snp2gene

	output:
		file 'snp2snp.tsv' into snp2snp_1, snp2snp_2, snp2snp_3

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)
	library(data.table)

	snp2gene <- read_tsv("${SNP2GENE}", col_types = "cc")
	read_tsv("${TAB2}", col_types = cols(.default = col_character())) %>%
		rename(gene_1 = "Official Symbol Interactor A", gene_2 = "Official Symbol Interactor B") %>%
		select(gene_1, gene_2) %>%
		data.table %>%
		merge(snp2gene, by.x = 'gene_1', by.y = 'gene', allow.cartesian = TRUE) %>%
		merge(snp2gene, by.x = 'gene_2', by.y = 'gene', allow.cartesian = TRUE, suffix = c('_1', '_2')) %>%
		mutate(uniq_snp_id = cbind(snp_1, snp_2) %>% apply(1, sort) %>% apply(2, paste, collapse = "_"),
			   uniq_gene_id = cbind(gene_1, gene_2) %>% apply(1, sort) %>% apply(2, paste, collapse = "_")) %>%
		write_tsv("snp2snp.tsv")
	"""

}

process preprocess_data {

	input:
		file BED from bed
		file bim
		file FAM from fam
		file SNP2SNP from snp2snp_1 
		file EXCLUDED from excluded_snps
		val I from params.nperm

	output:
		set 'post_qc.bed', 'post_qc.bim', 'post_qc.fam', 'pheno' into filtered_bed

	"""
	# QC + keep only SNPs in models + remove excluded SNPs
	cut -f3 ${SNP2SNP} >tmp
	cut -f4 ${SNP2SNP} >>tmp
	sort tmp | uniq >snps_in_models
	comm -3 snps_in_models <(sort ${EXCLUDED}) >included_snps
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
		file SNP2SNP from snp2snp_2

	output:
		set val(I), 'scored_interactions.regression.txt' into snp_pairs, snp_pairs_null 

	"""
	epistasis_regression.nf --bfile ${BED.baseName} --pheno ${PHENO} --i ${I} -profile cluster

	if [ `wc -l < scored_interactions.regression.txt` -gt "1" ]; then
		# exhaustive LD pruning
		cut -f1 scored_interactions.regression.txt >tmp
		cut -f2 scored_interactions.regression.txt >>tmp
		sort tmp | uniq >included_snps

		plink -bfile ${BED.baseName} -extract included_snps -allow-no-sex -r2
		sed 's/^ \\+//' plink.ld | sed 's/ \\+/\t/g' | sed 's/\t\$//' >ld.tsv

		R -e '
		library(tidyverse); 
		ld_ok <- read_tsv("ld.tsv", col_types = "idcidcd") %>%
		filter(R2 < 0.75) %>%
		mutate(uniq_snp_id = cbind(SNP_A, SNP_B) %>% apply(1, sort) %>% apply(2, paste, collapse = "_")) %>%
		.\$uniq_snp_id
		snp_models <- read_tsv("${SNP2SNP}") %>%
		.\$uniq_snp_id	
		read_tsv("scored_interactions.regression.txt", col_types = "ccccddd") %>%
		mutate(uniq_snp_id = cbind(SNP1, SNP2) %>% apply(1, sort) %>% apply(2, paste, collapse = "_")) %>%
		filter(uniq_snp_id %in% ld_ok & uniq_snp_id %in% snp_models) %>%
		write_tsv("scored_interactions.regression.txt")'
	else
		echo "`cat scored_interactions.regression.txt`\tuniq_snp_id" >header
		sed 's/\tP\t/\t0.0001\t/' header | sed 's/\tSTAT\t/\t0\t/' | sed 's/\tBETA_INT\t/\t0\t/' >vals
		cat header vals >scored_interactions.regression.txt
	fi
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
		file SNP2SNP from snp2snp_3 
		file TAB2 from tab2

	output:
		set val(I), 'scored_gene_pairs.tsv' into scored_gene_pairs, scored_gene_pairs_null
		file 'sign_snp_pairs.tsv' optional true into sign_snp_pairs

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)
        library(data.table)

	snp2snp <- read_tsv('$SNP2SNP', col_types = 'cccccc') %>%
		data.table
	snp2gene <- read_tsv('$SNP2GENE', col_types = 'cc')
	threshold <- lapply(list.files(pattern = 'permuted_association_'), function(x) {
			read_tsv(x, col_types = 'ccccddd') %>%
				mutate(uniq_snp_id = cbind(SNP1, SNP2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
				data.table %>%
                		merge(snp2snp, by = 'uniq_snp_id', allow.cartesian = TRUE) %>%
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

	if (${I} == 1) {
		separate(snp_pairs, uniq_snp_id, into = c('snp_1','snp_2'), sep = '_') %>%
			unique %>%
			write_tsv('sign_snp_pairs.tsv')
	}

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
				  tau_001 = prod(P[P < 0.001]),
			 	  snp_pairs = paste(unique(uniq_snp_id), collapse = ',')) %>%
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
		file SIGN_SNP_PAIRS from sign_snp_pairs
		file GENE_PAIRS from gene_pairs
		file TAB2 from tab2

	output:
		file SIGN_SNP_PAIRS
		file 'sign_gene_pairs.tsv'
		file 'sign_pathways.tsv'

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
		separate(uniq_gene_id, into = c('gene_1', 'gene_2'), sep = '_') %>%
		select(-experiment) %>%
		select(gene_1, gene_2, P, everything())
	write_tsv(gene_pairs_assoc, 'sign_gene_pairs.tsv')
		
	# compute pathway association
	net <- read_tsv('${TAB2}', col_types = cols(.default = col_character())) %>%
		select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
		graph_from_data_frame
 
	# TODO create background
	bg <- NULL

	lapply(filter(gene_pairs_assoc, P < 0.05)\$uniq_gene_id, function(uniq_gene_id) {
			genes <- unlist(strsplit('a_b', split = '_'))
			# TODO more than 1?
			delete_edges(gsub('_', '|', 'uniq_gene_id', fixed = TRUE)) %>%
				shortest_paths(from = genes[1], to = genes[2]) %>%
				goenrich(bg = bg)
	}) %>%
		do.call(bind_rows, .) %>%
		write_tsv('sign_pathways.tsv')	
	"""

}
