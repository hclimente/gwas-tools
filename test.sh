#!/bin/bash 
set -e

cd test
./test_impute.sh -with-docker hclimente/gwas-tools:latest
./test_qc_gwas.sh -with-docker hclimente/gwas-tools:latest
./test_snp2gene.sh -with-docker hclimente/gwas-tools:latest
./test_lift_coordinates.sh -with-docker hclimente/gwas-tools:latest
./test_run_vegas.sh -with-docker hclimente/gwas-tools:latest
