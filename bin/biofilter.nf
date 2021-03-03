#!/usr/bin/env nextflow

params.out = '.'

params.loki = ''
snp2gene = file(params.snp2gene)


if (params.loki == '') {

    process build_loki {

        publishDir "$params.out", overwrite: true, mode: "copy"

        output:
            file "loki_*.db" into loki

        """
        loki-build.py --verbose --knowledge loki_`date +%Y%m%d`.db --update
        """
    
    }

} else {
    loki = file(params.loki)
}

process get_genes {

    input:
        file SNP2GENE from snp2gene

    output:
        file 'genes' into genes

    """
    cut -f2 ${SNP2GENE} | tail -n +2 | sort | uniq >genes
    """

}

process get_gene_models {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file LOKI from loki
        file GENES from genes

    output:
        file 'biofilter.gene.models' into gene_models

    """
    cat << EOF >parameters.txt
    KNOWLEDGE           ${LOKI}
    GENE_FILE           ${GENES}
    MODEL               gene
    EOF

    biofilter.py parameters.txt
    """

}

gene_models. splitText( by: 1000, keepHeader: true, file: true ).set { gene_models_split }

process make_snp_models {

    input:
        file SNP_MODELS from gene_models_split
        file SNP2GENE from snp2gene

    output:
        file 'snp_models_split.tsv' into snp_models_split

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    snp2gene <- read_tsv("${SNP2GENE}")

    read_tsv("${SNP_MODELS}") %>%
        inner_join(snp2gene, by = c('#gene1' = 'gene')) %>%
        inner_join(snp2gene, by = c('gene2' = 'gene'), suffix = c("_1", "_2")) %>%
        filter(snp_1 != snp_2) %>%
        rename(gene_1 = `#gene1`, gene_2 = gene2) %>%
        select(snp_1, gene_1, snp_2, gene_2, `score(src-grp)`) %>%
        separate(`score(src-grp)`, into = c('num_sources', 'num_instances'), sep = '-') %>%
        write_tsv('snp_models_split.tsv', col_names = FALSE)
    """

}

process join_snp_models {

    publishDir "$params.out", overwrite: true, mode: "move"

    input:
        file 'snp_models_split*' from snp_models_split.collect()

    output:
        file 'snp_models.tsv'

    """
    echo "snp_1\tgene_1\tsnp_2\tgene_2\tnum_sources\tnum_instances" >snp_models.tsv
    cat snp_models_split* >>snp_models.tsv
    """
}
