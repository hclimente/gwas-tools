#!/usr/bin/env bash

wget --no-check-certificate https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0.tar.gz && R -e 'install.packages("dmGWAS_3.0.tar.gz", repos = NULL, type="source")'
wget --user-agent="Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1)" https://vegas2.qimrberghofer.edu.au/vegas2v2 && chmod a+x vegas2v2
R -e 'install.packages("igraph", repos = "http://cran.us.r-project.org")'
