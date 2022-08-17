#!/usr/bin/env nextflow

include { dmgwas_nf as dmgwas } from './dmgwas.nf'
include { heinz_nf as heinz } from './heinz.nf'
include { hotnet2_nf as hotnet2 } from './hotnet2.nf'
include { lean_nf as lean } from './lean.nf'
include { scones_old_nf as scones } from './scones.nf'
include { sigmod_nf as sigmod } from './sigmod.nf'
include { snp2gene_nf as snp2gene } from './snp2gene.nf'
include { snp_association_nf as snp_assoc } from './snp_association.nf'
include { vegas2_nf as gene_assoc } from './vegas2.nf'
include { get_bfile } from './templates/utils.nf'

// Parameters
/// splits
params.splits = 5

/// snp association
bfile = get_bfile(params.bfile)

/// vegas2
params.gencode = 28
params.genome = '38'
params.covars = ''
params.ld_controls = ''
params.buffer = 0

bfile_controls = get_bfile(params.bfile_controls)

/// edgelist selection
params.edgelist = null

//// dmgwas
params.dmgwas_r = 0.1
params.dmgwas_d = 2

params.heinz_fdr = 0.1

params.hotnet2_heat_permutations = 1000
params.hotnet2_lfdr_cutoff = 0.05
params.hotnet2_network_permutations = 100

params.scones_criterion = 'consistency'
params.scones_network = 'gi'
params.scones_score = 'chi2'

params.sigmod_lambdamax = 1
params.sigmod_nmax = 300
params.sigmod_maxjump = 10

process download_hint {

    output:
        file "hint_hgnc.tsv"

    script:
    template 'dbs/hint.sh'

}

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

process scones_genes {

    tag { SCONES_OUT.getBaseName() }

    input:
        path SCONES_OUT
        path SNP2GENE

    output:
        path "${SCONES_OUT.getBaseName()}.genes.tsv"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    snp2gene <- read_tsv("${SNP2GENE}")

    read_tsv("${SCONES_OUT}", col_types = "c") %>%
        inner_join(snp2gene, by = "snp") %>%
        select(gene) %>%
        unique %>%
        write_tsv("${SCONES_OUT.getBaseName()}.genes.tsv")
    """

}

process standardize_outputs {

    input:
        path OUT

    output:
        path "standardized.tsv"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(tools)

    ext <- file_ext("${OUT}")
    selected <- read_tsv("${OUT}")

    if ("Gene" %in% colnames(selected))
        selected <- rename(selected, gene = Gene)

    if ("PLEAN" %in% colnames(selected))
        selected <- filter(selected, PLEAN < 0.05)

    selected %>%
        mutate(method = sub(paste0('.', ext), '', '${OUT}'),
               method = sub('[^.]+.', '', method)) %>%
        select(gene, method) %>%
        write_tsv('standardized.tsv')
    """
}

process build_consensus {

    input:
        path "genes_*"

    output:
        path "stable_consensus.tsv"

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    list.files(pattern="genes_") %>%
        lapply(read_tsv) %>%
        bind_rows %>%
        group_by(gene, method) %>%
        summarize(n = n(),
                  tmp = paste0(method, "(", n(), ")")) %>%
        unique %>%
        group_by(gene) %>%
        summarize(n_selected = sum(n),
                  methods = paste(tmp, collapse = ',')) %>%
        arrange(-n_selected) %>%
        write_tsv("stable_consensus.tsv")
    """

}

workflow {
    main:

        if (params.edgelist == null) {
            download_hint()
            edgelist = download_hint.out
        } else {
            edgelist = file(params.edgelist)
        }

        split_data(bfile, 1..params.splits, params.splits)
        snp_assoc(split_data.out, params.covars)
        gene_assoc(snp_assoc.out, bfile_controls, params.gencode, params.genome, params.buffer, '')
        snp2gene(bfile[1], params.gencode, params.genome, params.buffer)
        dmgwas(gene_assoc.out, edgelist, params.dmgwas_d, params.dmgwas_r)
        heinz(gene_assoc.out, edgelist, params.heinz_fdr)
        hotnet2(gene_assoc.out, edgelist, params.hotnet2_lfdr_cutoff, params.hotnet2_network_permutations, params.hotnet2_heat_permutations)
        lean(gene_assoc.out, edgelist)
        scones(split_data.out, edgelist, params.scones_network, snp2gene.out, params.scones_score, params.scones_criterion, null, null)
        scones_genes(scones.out, snp2gene.out)
        sigmod(gene_assoc.out, edgelist, params.sigmod_lambdamax, params.sigmod_nmax, params.sigmod_maxjump)

        outputs = dmgwas.out
            .mix(heinz.out, hotnet2.out, lean.out, scones_genes.out, sigmod.out)
        standardize_outputs(outputs)
        build_consensus(standardize_outputs.out.collect())
    emit:
        build_consensus.out.collectFile(name: 'stable_consensus.tsv', storeDir: '.')
}
