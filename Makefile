PYTHON ?= python3
RSCRIPT ?= Rscript

.PHONY: all data models figures audit clean

all: data models figures audit

data:
	$(PYTHON) scripts/01_prepare_data.py

models:
	$(RSCRIPT) scripts/02_run_models_leakage_safe.R

figures:
	$(RSCRIPT) scripts/03_make_figures.R

audit:
	$(PYTHON) scripts/04_audit_project.py

clean:
	rm -f data/processed/*.csv
	rm -f results/tables/*.csv results/tables/*.json
	rm -f results/figures/*.png
