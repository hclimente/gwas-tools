FROM r-base

RUN apt-get update \
    && apt-get install -y wget unzip tar python sed bedtools libcurl4-gnutls-dev libxml2-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*
RUN mkdir /gwas-tools
ENV PATH="/gwas-tools:${PATH}"
WORKDIR /gwas-tools
RUN R -e "install.packages(c('mvtnorm', 'corpcor', 'tidyverse', 'magrittr'), repos = 'http://cran.us.r-project.org')"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('martini')"
RUN wget http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz \
    && tar -zxvf gtool_v0.7.5_x86_64.tgz
RUN wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz \
    && tar -zxvf impute_v2.3.2_x86_64_static.tgz && mv impute_v2.3.2_x86_64_static/impute2 .
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && chmod a+x liftOver
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && tar jxvf bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && cp bin/* .
RUN wget https://vegas2.qimrberghofer.edu.au/vegas2v2 && chmod a+x vegas2v2
RUN wget https://www.cog-genomics.org/static/bin/plink180913/plink_linux_x86_64.zip \
    && unzip plink_linux_x86_64.zip
RUN wget --no-check-certificate https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0.tar.gz \
    && R -e 'install.packages("dmGWAS_3.0.tar.gz", repos = NULL, type="source")'
WORKDIR /home/docker