#!/usr/bin/env bash

# GTOOL
wget http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz  && \
    tar -zxvf gtool_v0.7.5_x86_64.tgz && \
    rm -rf gtool_v0.7.5_x86_64.tgz LICENCE example/

# IMPUTE
wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
    tar -zxvf impute_v2.3.2_x86_64_static.tgz && \
    mv impute_v2.3.2_x86_64_static/impute2 . && \
    rm -rf impute_v2.3.2_x86_64_static* 

# MB-MDR
wget http://bio3.giga.ulg.ac.be/software/mbmdr-4.4.1/mbmdr-4.4.1-linux-64bits.out && \
    mv mbmdr-4.4.1-linux-64bits.out mbmdr && \
    chmod a+x mbmdr

# AntEpiSeeker
wget http://nce.ads.uga.edu/~romdhane/AntEpiSeeker/AntEpiSeeker1.0_linux.zip && \
    unzip AntEpiSeeker1.0_linux.zip && \
    mv AntEpiSeeker1.0_linux/AntEpiSeeker . && \
    rm -rf AntEpiSeeker1.0_linux*
