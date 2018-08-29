FROM ubuntu:18.04

RUN apt-get update \
    && apt-get install -y wget unzip tar python \
    && rm -rf /var/lib/apt/lists/*
RUN mkdir /gwas-tools
WORKDIR /gwas-tools
RUN wget https://www.cog-genomics.org/static/bin/plink180807/plink_linux_x86_64.zip \
    && unzip plink_linux_x86_64.zip
RUN wget http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz \
    && tar -zxvf gtool_v0.7.5_x86_64.tgz
RUN wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz \
    && tar -zxvf impute_v2.3.2_x86_64_static.tgz && mv impute_v2.3.2_x86_64_static/impute2 .
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && chmod a+x liftOver
ENV PATH="/gwas-tools:${PATH}"