#!/usr/bin/env nextflow

params.out = '.'

bim = file(params.bim)
params.gencode_version = 28
params.grch_version = '38'
params.buffer = 50000

process download_hgnc {

    output:
        file 'non_alt_loci_set.txt' into hgnc

    script:
    template 'dbs/hgnc.sh'
    
}

process download_annotation {

    input:
        val GENCODE_VERSION from params.gencode_version
        val GRCH_VERSION from params.grch_version

    output:
        file 'gff3' into gff
    
    script:
    template 'dbs/gencode.sh'

}

process extract_genes {

	input:
		file gff

	output:
		file 'genes.gff' into genes_gff

	"""
	awk '\$3 == "gene"' $gff >genes.gff
	"""

}

process gff2bed {

	input:
		file genes_gff
		val BUFFER from params.buffer

	output:
		file 'genes.bed' into genes_bed

	"""
	gff2bed < $genes_gff >tmp
	awk '{\$2 = \$2 - ${BUFFER}; \$2 = (\$2<0?0:\$2); \$3 = \$3 + ${BUFFER}; \$3 = (\$3<0?0:\$3); print}' OFS='\\t' tmp >genes.bed
	"""

}

process bim2bed {

	input:
		file BIM from bim

	output:
		file 'snps.bed' into snps_bed

	"""
	awk '{print "chr" \$1 "\\t" \$4 "\\t" \$4 "\\t" \$2 "\\t.\\t." }' ${BIM} >tmp
	sed 's/chr23/chrX/' tmp >snps.bed
	"""

}

process snp2gene {

	input:
		file snps_bed
		file genes_bed

	output:
		file 'snp2ensembl.tsv' into snp2ensembl

	"""
	bedtools intersect -a $snps_bed -b $genes_bed -wa -wb >tmp
	cut -f4,16 tmp | sed 's/\\tID=/\\t/' | sed 's/\\.[0-9]\\+;.\\+//' >snp2ensembl.tsv
	"""

}

process ensembl2hgnc {

	publishDir "$params.out", overwrite: true, mode: "copy"

	input:
		file hgnc
		file snp2ensembl

	output:
		file 'snp2hgnc.tsv' into snp2hgnc

	"""
	#!/usr/bin/env Rscript
	library(tidyverse)
	library(magrittr)
	ensembl2hgnc <- read_tsv('$hgnc') %>%
		select(symbol, ensembl_gene_id)
	read_tsv('$snp2ensembl', col_names = FALSE) %>%
		set_colnames(c('snp','ensembl_gene_id')) %>%
		inner_join(ensembl2hgnc, by = 'ensembl_gene_id') %>%
		select(snp, symbol) %>%
		rename(gene = symbol) %>%
		write_tsv('snp2hgnc.tsv')
	"""

}