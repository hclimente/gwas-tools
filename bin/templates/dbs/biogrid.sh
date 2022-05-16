# Getting the protein-protein interactions
wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab2.zip
unzip BIOGRID-MV-Physical-LATEST.tab2.zip
mv *.tab2.txt tab2
# if not at root of gwas-tools directory, change path so it gets into gwas-tools/test/data/
mv tab2 test/data/
rm BIOGRID-MV-Physical-LATEST.tab2.zip
