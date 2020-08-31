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
		set 'post_qc.bed', 'post_qc.bim', 'post_qc.fam', 'pheno' into filtered_bed_pipeline, filtered_bed_qc
		file 'snps_75.prune.in' into qc_snps

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

	tag { "${I}" }

	input:
		each I from 1..(params.nperm + 1)
		set file(BED), file(BIM), file(FAM), file(PHENO) from filtered_bed_pipeline
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
	
	tag { "${I}" }

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
				head(1) %>%
				select(P)
			}) %>%
		do.call(bind_rows, .) %>%
		.\$P %>%
		quantile(.05)

	snp_pairs <- read_tsv('${SNP_PAIRS}', col_types = 'ccccddd') %>%
		mutate(uniq_snp_id = cbind(SNP1, SNP2) %>% apply(1, sort) %>% apply(2, paste, collapse = '_')) %>%
		filter(P < threshold) %>%
		mutate(Padj = (P * .05) / threshold) %>%
		inner_join(snp2snp, by = 'uniq_snp_id') %>%
		select(uniq_snp_id, Padj)

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
		summarize(tau_05 = prod(Padj[Padj < 0.05]),
				  tau_01 = prod(Padj[Padj < 0.01]),
				  tau_001 = prod(Padj[Padj < 0.001]),
			 	  snp_pairs = paste(unique(uniq_snp_id), collapse = ',')) %>%
        write_tsv('scored_gene_pairs.tsv')
    """

}

// create two subsets
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
		file SNP2GENE from snp2gene
		file TAB2 from tab2
		file SNPS from qc_snps // TODO should be filtered bim?

	output:
		file 'sign_gene_pairs.tsv'
		file 'sign_pathways.tsv'

	"""
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(igraph)
    library(clusterProfiler)
    library(msigdbr)

	# compute gene-pair association
    gene_pairs <- read_tsv('${GENE_PAIRS}') %>%
        mutate(experiment = 'original')
    studied_gene_pairs <- gene_pairs\$uniq_gene_id		
    permutations <- list.files(pattern = 'permuted_association_')
    N <- length(permutations) + 1

    gene_pairs_assoc <- lapply(permutations, function(x) {
        read_tsv(x, col_types = 'cdddc') %>%
            filter(uniq_gene_id %in% studied_gene_pairs)
    }) %>%
        do.call(bind_rows, .) %>%
        mutate(experiment = 'permutation') %>%
        bind_rows(gene_pairs) %>%
        group_by(uniq_gene_id) %>%
        mutate(P_tau_05 = rank(tau_05) / N,
               P_tau_01 = rank(tau_01) / N, 
               P_tau_001 = rank(tau_001) / N,
               min_P = pmin(P_tau_05, P_tau_01, P_tau_001),
               P = rank(min_P) / N ) %>%
        filter(experiment == 'original') %>%
        separate(uniq_gene_id, into = c('gene_1', 'gene_2'), sep = '_') %>%
        select(-experiment) %>%
        select(gene_1, gene_2, P, everything())
    write_tsv(gene_pairs_assoc, 'sign_gene_pairs.tsv')
    
    sign_pairs <- filter(gene_pairs_assoc, P < .05)
    sign_genes <- c(sign_pairs\$gene_1, sign_pairs\$gene_2) %>% unique

    # compute pathway association
    net <- read_tsv('${TAB2}', col_types = cols(.default = col_character())) %>%
        select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
        graph_from_data_frame(directed = FALSE)
 
	# create background
    snps <- read_tsv('${SNPS}', col_names = FALSE)\$X1
    snp2gene <- read_tsv('${SNP2GENE}', col_types = 'cc')
    bg <- snp2gene\$gene[snp2gene\$snp %in% snps]

    # import MSigDB Gene Sets
    m_df = msigdbr(species = "Homo sapiens")
    m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

    mapply(function(gene_1, gene_2) {
        # TODO more than 1?
        delete_edges(net, paste(gene_1, gene_2, sep = '|')) %>%
            shortest_paths(from = gene_1, to = gene_2) %>%
            .\$vpath %>% unlist %>% names %>%
            intersect(sign_genes) %>%
            enricher(TERM2GENE = m_t2g, universe = bg, pAdjustMethod = 'bonferroni')
    }, sign_pairs\$gene_1, sign_pairs\$gene_2) %>%
        lapply(as_tibble) %>%
        bind_rows %>%
        write_tsv('sign_pathways.tsv')	
    """

}

if (params.prs != '') {
	process adjust_snp_pairs {

		publishDir "$params.out", overwrite: true, mode: "copy"

		input:
			file SIGN_SNP_PAIRS from sign_snp_pairs
			file PRS from file(params.prs)
			set file(BED), file(BIM), file(FAM), file(PHENO) from filtered_bed_qc

		output:
			file 'prs_adjusted_sign_snp_pairs.tsv'

		"""
		#!/usr/bin/env Rscript

		library(tidyverse)
		library(snpStats)

		snp_pairs <- read_tsv('${SIGN_SNP_PAIRS}', col_types = 'ccd')
		gwas <- read.plink("${BED}", "${BIM}", "${FAM}")
		prs <- read_tsv('${PRS}')\$prs

		X <- as(gwas[['genotypes']], 'numeric')
		X[X == 0] = 'AA'
		X[X == 1] = 'Aa'
		X[X == 2] = 'aa'
		y <- gwas[['fam']][['affected']]

		gwas <- as_tibble(X) %>%
			mutate(y = y, prs = prs)
		rm(X,y)

		prs_adjust <- function(snp_1, snp_2, ...) {
			df <- as.data.frame(gwas[, c('y', 'prs', snp_1, snp_2)])
			colnames(df) <- c('Y', 'PRS', 'SNP1', 'SNP2')
			PRS_adjusted = lm(Y ~ PRS + SNP1 + SNP2 + SNP1*SNP2, df, na.action=na.exclude)
        	        summary(PRS_adjusted)\$coefficients[5, 4]
		}

		snp_pairs %>%
			mutate(prs_adjusted_p =  pmap_dbl(snp_pairs, prs_adjust)) %>%
			write_tsv('prs_adjusted_sign_snp_pairs.tsv')
	"""
	}
}
