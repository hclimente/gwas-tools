###############################
# GLOBALS
CONDA_ENV = ./env/
CONDA_ACTIVATE = eval "$$(conda shell.bash hook)"; conda activate gwas-tools 
SHELL=bash

.PHONY: $(CONDA_ENV) clean conda docker setup

###############################
# COMMANDS
setup: $(CONDA_ENV)
	$(CONDA_ACTIVATE); mamba install --force python-lsp-server

$(CONDA_ENV): environment.yml
	mamba env create --force --prefix $(CONDA_ENV) --file environment.yml

conda: environment.yml
	mamba env create --force --file environment.yml
	$(CONDA_ACTIVATE); bash nonconda_deps.sh
	mv vegas2v2 bin/

docker:
	docker build -t hclimente/gwas-tools .
	docker build -t hclimente/gwas-tools-extra - <Dockerfile_extra

clean:
	rm -rf env/

