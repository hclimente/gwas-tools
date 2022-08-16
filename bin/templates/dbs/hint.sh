wget http://hint.yulab.org/download/HomoSapiens/htb/hq/ -O hint.tsv
cut -f3,4 hint.tsv | sed 's/Gene_A/gene1/' | sed 's/Gene_B/gene2/' >hint_hgnc.tsv
