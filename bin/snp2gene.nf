#!/usr/bin/env nextflow

params.out = '.'

params.gencode_version = 41
params.genome_version = '38'
params.buffer = 50000

process download_hgnc {

    output:
        path 'non_alt_loci_set.txt'

    script:
    template 'dbs/hgnc.sh'
    
}

process download_gencode {

    input:
        val GENCODE_VERSION
        val GRCH_VERSION

    output:
        path 'gff3'
    
    script:
    template 'dbs/gencode.sh'

}

process gff2bed {

	input:
		path GFF
		val BUFFER

	output:
		path 'genes.bed'

	"""
	awk '\$3 == "gene"' ${GFF} >genes.gff
	gff2bed <genes.gff >tmp
	awk '{\$2 = \$2 - ${BUFFER}; \$2 = (\$2<0?0:\$2); \$3 = \$3 + ${BUFFER}; \$3 = (\$3<0?0:\$3); print}' OFS='\\t' tmp >genes.bed
	"""

}

process bim2bed {

	input:
		path BIM

	output:
		path 'snps.bed'

	"""
	awk '{print "chr" \$1 "\\t" \$4 "\\t" \$4 "\\t" \$2 "\\t.\\t." }' ${BIM} >tmp
	sed 's/chr23/chrX/' tmp >snps.bed
	"""

}

process map_snps_to_genes {

	input:
		path snps_bed
		path genes_bed

	output:
		path 'snp2ensembl.tsv'

	"""
	bedtools intersect -a $snps_bed -b $genes_bed -wa -wb >tmp
	cut -f4,16 tmp | sed 's/\\tID=/\\t/' | sed 's/\\.[0-9]\\+;.\\+//' >snp2ensembl.tsv
	"""

}

process ensembl2hgnc {

    publishDir params.out, mode: 'copy'

	input:
		path snp2ensembl
		path hgnc

	output:
		path 'snp2hgnc.tsv'

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

workflow snp2gene_nf {
    take:
        bim
        gencode_version
        genome_version
        buffer
    main:
        download_hgnc()
        download_gencode(gencode_version, genome_version)
        gff2bed(download_gencode.out, buffer)
        bim2bed(bim)
        map_snps_to_genes(bim2bed.out, gff2bed.out)
        ensembl2hgnc(map_snps_to_genes.out, download_hgnc.out)
    emit:
        ensembl2hgnc.out
}

workflow {
    main:
        snp2gene_nf(file(params.bim), params.gencode_version, params.genome_version, params.buffer)
    emit:
        snp2gene_nf.out
}
