FROM r-base

RUN apt-get update \
    && apt-get install -y wget unzip tar python2 python3 python3-pip  sed bedtools libcurl4-gnutls-dev libxml2-dev libssl-dev gawk \
    && rm -rf /var/lib/apt/lists/*
RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2 get-pip.py \
    && rm get-pip.py
RUN mkdir /gwas-tools
ENV PATH="/gwas-tools:${PATH}"
WORKDIR /gwas-tools
RUN R -e "install.packages(c('mvtnorm', 'corpcor', 'tidyverse', 'magrittr', 'BiocManager', 'cowplot'), repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install(c('martini','BioNet','twilight'))"
#RUN wget http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz \
#    && tar -zxvf gtool_v0.7.5_x86_64.tgz
#RUN wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz \
#    && tar -zxvf impute_v2.3.2_x86_64_static.tgz && mv impute_v2.3.2_x86_64_static/impute2 .
#RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && chmod a+x liftOver
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && tar jxvf bedops_linux_x86_64-v2.4.35.tar.bz2 \
    && cp bin/* .
RUN wget https://github.com/YuanlongLiu/SigMod/blob/20c561876d87a0faca632a6b93882fcffd719b17/SigMod_v2.zip \
    && unzip SigMod_v2.zip
RUN wget https://github.com/raphael-group/hotnet2/archive/refs/tags/v1.2.1.tar.gz \
    && tar xzf v1.2.1.tar.gz \
    && rm v1.2.1.tar.gz \
    && mv hotnet2-1.2.1 hotnet2 \
    && pip2 install -r hotnet2/requirements.txt
RUN wget --no-check-certificate https://bioinfo.uth.edu/dmGWAS/dmGWAS_3.0.tar.gz \
    && R -e 'install.packages("dmGWAS_3.0.tar.gz", repos = NULL, type="source")'
#RUN wget http://bio3.giga.ulg.ac.be/software/mbmdr-4.4.1/mbmdr-4.4.1-linux-64bits.out \
#    && mv mbmdr-4.4.1-linux-64bits.out mbmdr && chmod a+x mbmdr
#RUN wget --no-check-certificate https://ritchielab.org/files/RL_software/biofilter-2.4.1.tar.gz \
#    && tar -xvzf biofilter-2.4.1.tar.gz \
#    && pip install apsw \
#    && cd biofilter-2.4.1 && python2 setup.py install && cd ..
#RUN wget https://sites.fas.harvard.edu/\~junliu/BEAM/BEAM_linux.tar \
#    && tar -xvf BEAM_linux.tar && apt-get update && apt-get install -y lib32z1 lib32stdc++6
#RUN wget http://nce.ads.uga.edu/~romdhane/AntEpiSeeker/AntEpiSeeker1.0_linux.zip \
#    && unzip AntEpiSeeker1.0_linux.zip && mv AntEpiSeeker1.0_linux/AntEpiSeeker .
#RUN apt-get update && apt-get install -y procps
# Doesn't exist anymore - seems to be usefull for epistasis analysis only
#RUN wget https://gwas.biosciencedbc.jp/SNPInterForest/rf && chmod a+x rf && mv rf SNPInterForest
RUN R -e "install.packages(c('ranger','SKAT','biglasso','bigmemory','igraph','LEANR','CASMAP', 'doMC'), repos = 'http://cran.us.r-project.org')"
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip \
    && unzip plink_linux_x86_64_20200616.zip
#RUN wget https://vegas2.qimrberghofer.edu.au/vegas2v2
RUN wget --user-agent="Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1; SV1)" https://vegas2.qimrberghofer.edu.au/vegas2v2 \ 
    && chmod a+x vegas2v2
WORKDIR /home/docker
