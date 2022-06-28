#!/usr/bin/env nextflow

include { heinz_workflow } from './heinz.nf'
include { vegas2_workflow } from './vegas2.nf'

params.out = '.'
params.splits = 5

// gwas
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

print(bed.baseName)
algorithms = [
              'dmgwas',
              'heinz',
            //  'hotnet2',
              'lean',
            //  'scones',
              'sigmod'
            ]

// annotation
tab2 = file(params.tab2)
snp2gene = file(params.snp2gene)

// code
hotnet2_path = file('../hotnet2/hotnet2')

process split_data {

    input:
        file FAM
        each I 
        val K

    output:
        file 'samples_*.txt'

    script:
    template 'genotypes/train_test_split.sh'

}

process vegas {

    input:
        tuple path(BED), path(FAM), path(BIM)
        each SAMPLES

    output:
        file 'scores.txt'

    """
    plink --bfile ${BED.baseName} --keep ${SAMPLES} --make-bed --out input
    vegas2.nf --bfile input --genome 37 --gencode 31 --buffer 50000 --vegas_params '-top 10' -profile bigmem
    
    R --random-flags <<RSCRIPT
    library(tidyverse);
    download.file("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt", "ensg2hgnc.tsv")
    ensembl2hgnc <- read_tsv("ensg2hgnc.tsv") %>%
        select(symbol, ensembl_gene_id)
    vegas <- read_tsv(scored_genes.vegas.txt], col_types = 'iciddddddcd') %>%
        inner_join(ensembl2hgnc, by = c('Gene' = 'ensembl_gene_id')) %>%
        mutate(Gene = symbol,
               Pvalue = `Top-0.1-pvalue`) %>%
        select(Gene, Pvalue) %>%
        write_tsv(scores.txt)
    RSCRIPT
    """

}

process network_selection {
    
    input:
        each METHOD
        file SCORES
        file TAB2
        
    output:
        file "selected.txt"
        
    """
    ${METHOD}.nf --scores ${SCORES} --tab2 ${TAB2}
    """

}

process build_consensus {

    input:
        file "selected_*.txt"

    output:
        file "stable_consensus.txt"

    """
    """

}

workflow {
    main:
        split_data(fam, 1..params.splits, params.splits)
        vegas(tuple(bed, fam, bim), split_data.out)
        // heinz_workflow(vegas.out, tab2, 0.1)
        // build_consensus(network_selection.out.collectFile)
    emit:
        split_data.out
}
