#!/usr/bin/env nextflow

include { download_hgnc; download_gencode } from './snp2gene.nf'
include { get_bfile } from './templates/utils.nf'

params.buffer = 0
params.gencode = 28
params.genome = '38'
params.snp_association = ''
params.vegas2_params = ''

bfile_ld_controls = get_bfile(params.bfile_ld_controls)

process compute_gene_coords {

    input:
        path GFF
        val BUFF

    output:
        path 'glist_ensembl'

    """
    awk '\$3 == "gene"' $GFF >genes.gff
    gff2bed < genes.gff | cut -f1-4 | sed 's/\\.[^\\t]\\+\$//' | sed 's/^chr//' >tmp
    awk '{\$2 = \$2 - ${BUFF}; \$3 = \$3 + ${BUFF}} 1' tmp | awk '\$2 < 0 {\$2 = 0} 1' >buffered_genes
    sed 's/^XY/25/' buffered_genes | sed 's/^X/23/' | sed 's/^Y/24/' | sed 's/^M/26/' | awk '\$1 <= 24' >glist_ensembl
    """

}

process vegas2 {

    tag { "${SNP_ASSOC}, chr ${CHR}" }

    input:
        path SNP_ASSOC
        tuple path(BED), path(BIM), path(FAM)
        path GLIST
        val VEGAS_PARAMS
        each CHR

    output:
        tuple val("${SNP_ASSOC.getBaseName()}"), path("chr_${CHR}.tsv")

    """
    cut -f3,4 ${SNP_ASSOC} | tail -n +2 >assoc

    vegas2v2 -G -snpandp assoc -custom `pwd`/${BED.baseName} -glist ${GLIST} -out scored_genes -chr $CHR ${VEGAS_PARAMS}
    sed 's/"//g' scored_genes.out | sed 's/ /\\t/g' >tmp
    R -e 'library(tidyverse); read_tsv("tmp", col_types = "iciddddddcd") %>% filter(!duplicated(Gene)) %>% write_tsv("chr_${CHR}.tsv")'
    """

}

process merge_chromosomes {

    tag { KEY }

    input:
        tuple val(KEY), path('chr_*')
		path HGNC

    output:
        path "${KEY}.vegas2.tsv"

    """
	#!/usr/bin/env Rscript

    library(tidyverse)

	ensembl2hgnc <- read_tsv('$HGNC') %>%
		select(symbol, ensembl_gene_id)

    list.files(pattern="chr_") %>%
        lapply(read_tsv, col_types = 'cciiiiddcd') %>%
        bind_rows %>%
        inner_join(ensembl2hgnc, by = c("Gene" = "ensembl_gene_id")) %>%
        mutate(Gene = symbol) %>%
        select(-symbol) %>%
        write_tsv('${KEY}.vegas2.tsv')
    """

}

workflow vegas2_nf {
    take:
        snp_association
        bfile_ld_controls
        gencode_version
        genome_version
        buffer
        vegas2_params
    main:
        download_gencode(gencode_version, genome_version)
        compute_gene_coords(download_gencode.out, params.buffer)

        chromosomes = snp_association
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .unique()

        vegas2(snp_association, bfile_ld_controls, compute_gene_coords.out, vegas2_params, chromosomes)

        vegas_outputs = vegas2.out
            .groupTuple()

        download_hgnc()
        merge_chromosomes(vegas_outputs, download_hgnc.out)
    emit:
        merge_chromosomes.out
}

workflow {
    main:
        vegas2_nf(file(params.snp_association), bfile_ld_controls, params.gencode, params.genome, params.buffer, params.vegas2_params)
    emit:
        vegas2_nf.out
}
