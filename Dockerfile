FROM continuumio/miniconda3

# setup conda virtual environment
RUN conda install mamba -n base -c conda-forge
COPY environment.yml .
RUN mamba env update -n base -f environment.yml

# additional software
RUN mkdir /tools
ENV PATH="/tools:${PATH}"
WORKDIR /tools

COPY extra_deps.sh .
RUN bash extra_deps.sh && rm extra_deps.sh

# set up workdir
RUN mkdir /work
WORKDIR /work
