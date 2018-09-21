# Input variables:
#    - GENCODE_VERSION
#    - GENOME_VERSION
# Output file:
#    - gff3

if [ '$GENOME_VERSION' == 'GRCh37' ];
then
    tag1='/GRCh37_mapping'
    tag2='lift37'
else
    tag1=''
    tag2=''
fi

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}\$tag1/gencode.v${GENCODE_VERSION}\${tag2}.annotation.gff3.gz
gunzip -c *gff3.gz >gff3