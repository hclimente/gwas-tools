#!/usr/bin/env nextflow

params.out = '.'

process lean {

    publishDir params.out, mode: 'copy'
    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path EDGELIST

    output:
        path "${SCORES.getBaseName()}.lean.tsv"

    script:
    """
#!/usr/bin/env Rscript
library(igraph)
library(LEANR)
library(tidyverse)

# read edgelist
net <- read_tsv("${EDGELIST}") %>%
    graph_from_data_frame(directed = FALSE)

# read gene scores
gene_scores <- read_tsv('$SCORES') %>%
    filter(gene %in% names(V(net)))

scores <- gene_scores[['pvalue']]
names(scores) <- gene_scores[['gene']]

# run lean
results <- run.lean(scores, net, n_reps = 10000, 
                    add.scored.genes = TRUE, verbose = TRUE)

as_tibble(results\$restab) %>%
    mutate(gene = rownames(results\$restab)) %>%
    select(gene, everything()) %>%
    write_tsv('${SCORES.getBaseName()}.lean.tsv')
    """

}

workflow lean_nf {
    take:
        scores
        edgelist
    main:
        lean(scores, edgelist)
    emit:
        lean.out
}

workflow {
    main:
        lean_nf(file(params.scores), file(params.edgelist))
    emit:
        lean_nf.out
}
