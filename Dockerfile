FROM continuumio/miniconda3

# setup conda virtual environment
RUN conda install mamba -n base -c conda-forge
COPY environment.yml .
RUN mamba env update -n base -f environment.yml

# additional software
RUN mkdir /tools
ENV PATH="/tools:${PATH}"
WORKDIR /tools

COPY nonconda_deps.sh .
RUN bash nonconda_deps.sh && rm nonconda_deps.sh

# set up workdir
RUN mkdir /work
WORKDIR /work
