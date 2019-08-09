#!/usr/bin/env nextflow

params.out = '.'

params.gencode = 28
params.genome = '38'
params.vegas_params = ''
params.covar = ''
params.snp_association = ''
params.ld_controls = ''
params.buffer = 0

//////////////////////////////////////////////
///RUN VEGAS         ///
//////////////////////////////////////////////
process download_gencode {

    input:
        val GENCODE_VERSION from params.gencode
        val GRCH_VERSION from params.genome

    output:
        file 'gff3' into gff
    
    script:
    template 'dbs/gencode.sh'

}

process make_glist {

    input:
        file gff
        val BUFF from params.buffer

    output:
        file 'glist_ensembl' into glist

    """
    awk '\$3 == "gene"' $gff >genes.gff
    gff2bed < genes.gff | cut -f1-4 | sed 's/\\.[^\\t]\\+\$//' | sed 's/^chr//' >tmp
    awk '{\$2 = \$2 - ${BUFF}; \$3 = \$3 + ${BUFF}} 1' tmp | awk '\$2 < 0 {\$2 = 0} 1' >buffered_genes
    sed 's/^XY/25/' buffered_genes | sed 's/^X/23/' | sed 's/^Y/24/' | sed 's/^M/26/' | awk '\$1 <= 24' >glist_ensembl
    """

}

//////////////////////////////////////////////
///          SNP ASSOCIATION TEST          ///
//////////////////////////////////////////////
if (params.snp_association == ''){

  bed = file("${params.bfile}.bed")
  bim = file("${bed.getParent()}/${bed.getBaseName()}.bim")
  fam = file("${bed.getParent()}/${bed.getBaseName()}.fam")

  if (params.covar == '') {

      process compute_chisq {

          input:
              file BED from bed
              file BIM from bim
              file FAM from fam

          output:
              file 'snp_association' into snp_association

          """
          plink --bed ${BED} --bim ${BIM} --fam ${FAM} --assoc 
          awk 'NR > 1 && \$9 != "NA" { print \$2,\$9 }' OFS='\\t' plink.assoc  >snp_association
          """

      }

  } else {

    covar = file(params.covar)

    process regress_phenotypes_with_covars {

        input:
            file BED from bed
            file BIM from bim
            file FAM from fam
            file COVAR from covar

        output:
            file 'snp_association' into snp_association

        """
        plink --bed ${BED} --bim ${BIM} --fam ${FAM} --logistic --covar ${COVAR}
        awk 'NR > 1 && \$5 == "ADD" && \$9 != "NA" { print \$2,\$9 }' OFS='\\t' plink.assoc.logistic >snp_association
        """

    }
  }
} else {

  snp_association = file(params.snp_association)

}

//////////////////////////////////////////////
///         EXTRACT CONTROLS FOR LD        ///
//////////////////////////////////////////////
if (params.ld_controls == '') {
  process extract_controls {

      input:
          file BED from bed
          file bim
          file fam

      output:
          file 'plink.bed' into bed_controls
          file 'plink.bim' into bim_controls
          file 'plink.fam' into fam_controls

      """
      plink --bfile ${BED.baseName} --filter-controls --make-bed
      """

  }
} else {

  bed_controls = file("${params.ld_controls}.bed")
  bim_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.bim")
  fam_controls = file("${bed_controls.getParent()}/${bed_controls.getBaseName()}.fam")

}

//////////////////////////////////////////////
///                RUN VEGAS               ///
//////////////////////////////////////////////
process vegas {

    input:
        file BED from bed_controls
        file bim_controls
        file fam_controls
        file SNPASSOCIATION from snp_association
        file GLIST from glist 
        val VEGAS_PARAMS from params.vegas_params

    output:
        file 'scored_genes.vegas.duplicates.txt' into vegas_out

    """
    vegas2v2 -G -snpandp ${SNPASSOCIATION} -custom `pwd`/${BED.baseName} -glist ${GLIST} -out scored_genes ${VEGAS_PARAMS}
    sed 's/"//g' scored_genes.out | sed 's/ /\\t/g' >tmp
    mv tmp scored_genes.vegas.duplicates.txt
    """

}

process fix_vegas_output {

    publishDir "$params.out", overwrite: true, mode: "copy"

    input:
        file VEGAS from vegas_out

    output:
        file 'scored_genes.vegas.txt'

    """
    #!/usr/bin/env Rscript
    library(tidyverse)

    read_tsv('$VEGAS', col_types = 'iciddddddcd') %>%
        filter(!duplicated(Gene)) %>%
        write_tsv('scored_genes.vegas.txt')
    """

}
