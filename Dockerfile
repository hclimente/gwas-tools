FROM r-base

RUN apt-get update \
    && apt-get install -y wget unzip tar python python3-pip python2 sed bedtools libcurl4-gnutls-dev libxml2-dev libssl-dev gawk \
    && rm -rf /var/lib/apt/lists/*
RUN mkdir /gwas-tools
ENV PATH="/gwas-tools:${PATH}"
WORKDIR /gwas-tools
RUN R -e "install.packages(c('mvtnorm', 'corpcor', 'tidyverse', 'magrittr', 'BiocManager'), repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install(c('martini','BioNet'))"
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && tar jxvf bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && cp bin/* .
RUN wget --no-check-certificate https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0.tar.gz \
    && R -e 'install.packages("dmGWAS_3.0.tar.gz", repos = NULL, type="source")'
RUN apt-get update && apt-get install -y procps
RUN R -e "install.packages(c('ranger','SKAT','biglasso','bigmemory','igraph','LEANR','CASMAP'), repos = 'http://cran.us.r-project.org')"
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200219.zip \
    && unzip plink_linux_x86_20200219.zip
WORKDIR /home/docker
