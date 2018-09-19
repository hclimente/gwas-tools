# Input variables:
#    - GENCODE_VERSION
#    - TAG1
#    - TAG2
# Output file:
#    - gff3
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}${TAG1}/gencode.v${GENCODE_VERSION}${TAG2}.annotation.gff3.gz
gunzip -c *gff3.gz >gff3