#!/usr/bin/env nextflow

include { dmgwas_nf as dmgwas } from './dmgwas.nf'
include { heinz_nf as heinz } from './heinz.nf'
include { lean_nf as lean } from './lean.nf'
include { scones_nf as scones } from './scones.nf'
include { sigmod_nf as sigmod } from './sigmod.nf'
include { snp_association_nf as snp_assoc } from './snp_association.nf'
include { vegas2_nf as gene_assoc } from './vegas2.nf'
include { get_bfile } from './templates/utils.nf'

// Parameters
//////////////////////////////////////////////
// splits
params.splits = 5

// snp association
bfile = get_bfile(params.bfile)

// vegas2
params.gencode = 28
params.genome = '38'
params.covars = ''
params.ld_controls = ''
params.buffer = 0

bfile_controls = get_bfile(params.bfile_controls)

// network selection
tab2 = file(params.tab2)
snp2gene = file(params.snp2gene)

// dmgwas
params.r = 0.1
params.d = 2

// heinz
params.fdr = 0.1

// hotnet2
hotnet2_path = file('../hotnet2/hotnet2')

// scones
params.scones_criterion = 'consistency'
params.scones_network = 'gs'
params.scones_score = 'chi2'

// sigmod
// TODO remove
/* params.lambdamax = 1 */
/* params.nmax = 300 */
/* params.maxjump = 10 */
params.lambdamax = 2
params.nmax = 1
params.maxjump = 1

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
        /* dmgwas(gene_assoc.out, tab2, params.d, params.r) */
        heinz(gene_assoc.out, tab2, params.fdr)
        lean(gene_assoc.out, tab2)
        /* scones(bfile, tab2, params.scones_network, params.snp2gene, params.snp2gene, params.scones_score, params.scones_criterion) */
        /* map_scones_to_gene(scones.out) */
        sigmod(gene_assoc.out, tab2, params.lambdamax, params.nmax, params.maxjump)

        /* outputs = dmgwas.out */
        /*     .mix(heinz.out, lean.out, sigmod.out) */
        /*     .collect() */
        outputs = heinz.out
            .mix(lean.out, sigmod.out)
            .collect()
        build_consensus(outputs)
    emit:
        split_data.out
}
