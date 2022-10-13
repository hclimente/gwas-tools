#!/usr/bin/env nextflow

params.out = '.'

// sigmod params
params.lambdamax = 1
params.nmax = 300
params.maxjump = 10

process download_sigmod {

    output:
        path "SigMod_v2"

    """
    wget https://github.com/YuanlongLiu/SigMod/raw/20c561876d87a0faca632a6b93882fcffd719b17/SigMod_v2.zip && unzip SigMod_v2.zip
    """
}

process sigmod {

    publishDir params.out, mode: 'copy'
    tag { SCORES.getBaseName() }

    input:
        path SCORES
        path EDGELIST
        val LAMBDAMAX
        val NMAX
        val MAXJUMP
        path SIGMOD_PATH

    output:
        path "${SCORES.getBaseName()}.sigmod.txt"

    script:
    """
#!/usr/bin/env Rscript

require(igraph)
library(tidyverse)

scripts <- list.files('${SIGMOD_PATH}/R', pattern='*.R\$', full.names=TRUE, ignore.case=TRUE)
sapply(scripts, source, .GlobalEnv)

# read network
net <- read_tsv("${EDGELIST}") %>%
    filter(gene1 != gene2) %>%
    as.data.frame

# read vegas output
scores <- read_tsv('${SCORES}') %>% 
    rename(p = pvalue) %>%
    select(gene, p) %>%
    as.data.frame

# check weight_index = NULL
scored_net <- construct_scored_net(net, interaction_indices = c(1,2), gene_ps = scores)
res_info <- SigMod_bisection(net = scored_net, lambda_max = ${LAMBDAMAX}, nmax = ${NMAX}, maxjump = ${MAXJUMP})

save(scored_net, res_info, file = 'sigmod.RData')

data.frame(gene = names(V(res_info\$opt_module[[1]]))) %>%
    write_tsv("${SCORES.getBaseName()}.sigmod.txt")
    """

}

workflow sigmod_nf {
    take:
        scores
        edgelist
        lambdamax
        nmax
        maxjump
    main:
        download_sigmod()
        sigmod(scores, edgelist, lambdamax, nmax, maxjump, download_sigmod.out)
    emit:
        sigmod.out
}

workflow {
    main:
        sigmod_nf(file(params.scores), file(params.edgelist), params.lambdamax, params.nmax, params.maxjump)
    emit:
        sigmod_nf.out
}
