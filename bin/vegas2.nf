#!/usr/bin/env nextflow

params.gencode = 28
params.genome = '38'
params.vegas_params = ''
params.covar = ''
params.snp_association = ''
params.ld_controls = ''
params.buffer = 0
nextflow.enable.dsl=2

bed_controls = file("${params.ld_controls}.bed")
bim_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.bim")
fam_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.fam")
bfile_controls = tuple(bed_controls, bim_controls, fam_controls)

/// Create GLIST
//////////////////////////////////////////////
process download_gencode {

    input:
        val GENCODE_VERSION
        val GRCH_VERSION

    output:
        file 'gff3'
    
    script:
    template 'dbs/gencode.sh'

}

process compute_gene_coords {

    input:
        file GFF
        val BUFF

    output:
        file 'glist_ensembl'

    """
    awk '\$3 == "gene"' $GFF >genes.gff
    gff2bed < genes.gff | cut -f1-4 | sed 's/\\.[^\\t]\\+\$//' | sed 's/^chr//' >tmp
    awk '{\$2 = \$2 - ${BUFF}; \$3 = \$3 + ${BUFF}} 1' tmp | awk '\$2 < 0 {\$2 = 0} 1' >buffered_genes
    sed 's/^XY/25/' buffered_genes | sed 's/^X/23/' | sed 's/^Y/24/' | sed 's/^M/26/' | awk '\$1 <= 24' >glist_ensembl
    """

}

/// Vegas2
//////////////////////////////////////////////
process vegas {

    tag { "chr ${CHR}" }

    input:
        tuple file(BED), file(BIM), file(FAM)
        file SNPASSOCIATION
        file GLIST
        val VEGAS_PARAMS
        each CHR

    output:
        file 'scored_genes.vegas.txt'

    """
    vegas2v2 -G -snpandp ${SNPASSOCIATION} -custom `pwd`/${BED.baseName} -glist ${GLIST} -out scored_genes -chr $CHR ${VEGAS_PARAMS}
    sed 's/"//g' scored_genes.out | sed 's/ /\\t/g' >tmp
    R -e 'library(tidyverse); read_tsv("tmp", col_types = "iciddddddcd") %>% filter(!duplicated(Gene)) %>% write_tsv("scored_genes.vegas.txt")'
    """

}

process merge_chromosomes {

    input:
        file 'scored_genes_chr*'

    output:
        file 'scored_genes.vegas.txt'

    """
    head -n1 scored_genes_chr1 >scored_genes.vegas.txt
    tail -n +2 -q scored_genes_chr* | sort -n >>scored_genes.vegas.txt
    """

}

workflow vegas2_nf {
    take:
        snp_association
        bfile_controls
        gencode_version
        genome_version
        buffer
        vegas_params
    main:
        download_gencode(gencode_version, genome_version)
        compute_gene_coords(download_gencode.out, params.buffer)

        chromosomes = compute_gene_coords.out
            .splitCsv(header: false, sep: ' ')
            .map { row -> row[0] }
            .unique()

        vegas(bfile_controls, snp_association, compute_gene_coords.out, vegas_params, chromosomes)
        merge_chromosomes(vegas.out.collect())
    emit:
        merge_chromosomes.out
}

workflow {
    main:
        vegas2_nf(params.snp_association, bfile_controls, params.gencode, params.genome, params.buffer, params.vegas_params)
    emit:
        vegas2_nf.out
}
