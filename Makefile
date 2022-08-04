###############################
# GLOBALS
CONDA_ENV = ./env/
CONDA_ACTIVATE = eval "$$(conda shell.bash hook)"; conda activate $(CONDA_ENV)
SHELL=bash

.PHONY: $(CONDA_ENV) clean docker setup

###############################
# COMMANDS
setup: $(CONDA_ENV)
	$(CONDA_ACTIVATE); mamba install python-lsp-server

$(CONDA_ENV): environment.yml
	mamba env create --force --prefix $(CONDA_ENV) --file environment.yml

docker:
	docker build -t hclimente/gwas-tools .
	docker build -t hclimente/gwas-tools-extra - <Dockerfile_extra

clean:
	rm -rf env/

