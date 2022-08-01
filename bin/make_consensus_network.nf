#!/usr/bin/env nextflow

include { heinz_nf as heinz } from './heinz.nf'
include { lean_nf as lean } from './lean.nf'
include { sigmod_nf as sigmod } from './sigmod.nf'
include { snp_association_nf as snp_assoc } from './snp_association.nf'
include { vegas2_nf as gene_assoc } from './vegas2.nf'

// Parameters
//////////////////////////////////////////////

// splits
params.splits = 5

// snp association
bed = file("${params.bfile}.bed")
bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")
bfile = tuple(bed, bim, fam)

// vegas2
params.gencode = 28
params.genome = '38'
params.covars = ''
params.ld_controls = ''
params.buffer = 0

bed_controls = file("${params.bfile_controls}.bed")
bim_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.bim")
fam_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.fam")
bfile_controls = tuple(bed_controls, bim_controls, fam_controls)

// network selection
tab2 = file(params.tab2)
snp2gene = file(params.snp2gene)

// heinz
params.fdr = 0.1

// hotnet2
hotnet2_path = file('../hotnet2/hotnet2')

// sigmod
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10

// Pipeline
//////////////////////////////////////////////
process split_data {

    input:
        tuple path(BED), path(BIM), path(FAM)
        each I 
        val K

    output:
        tuple path("train_${I}.bed"), path("train_${I}.bim"), path("train_${I}.fam")

    script:
    template 'genotypes/train_test_split.sh'

}

process build_consensus {

    input:
        file "selected_*.txt"

    output:
        file "stable_consensus.tsv"

    """
#!/usr/bin/env Rscript

library(tidyverse)

list.files(pattern="selected_") %>%
    lapply(read_tsv) %>%
    bind_rows %>%
    group_by(gene) %>%
    summarize(n = n()) %>%
    write_tsv("stable_consensus.tsv")
    """

}

workflow {
    main:
        split_data(bfile, 1..params.splits, params.splits)
        snp_assoc(split_data.out, params.covars)
        gene_assoc(snp_assoc.out, bfile_controls, params.gencode, params.genome, params.buffer, '')
        heinz(gene_assoc.out, tab2, params.fdr)
        lean(gene_assoc.out, tab2)
        // sigmod(gene_assoc.out, tab2, params.lambdamax, params.nmax, params.maxjump)

        outputs = heinz.out
            .mix(lean.out)
            .collect()
        build_consensus(outputs)
    emit:
        split_data.out
}
