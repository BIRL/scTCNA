# Developer shortcuts for scTCNA. Everything runs inside the project virtualenv.
#
#   make venv install      create .venv and install the pinned environment
#   make data              download + convert every enabled dataset (data/raw)
#   make pipeline          process datasets into data/processed
#   make aggregate         build the app tables in data/app
#   make app               run the web app locally with auto-reload
#   make test / test-e2e   unit tests / Playwright browser tests
#   make lint / format     ruff check + format --check / apply formatting
#   make docker-build/run  build and run the deployment image locally
#
# Overridable: APP_DATA=data/app_syn (serve the synthetic tables), DATASETS="LUAD_Kim2020 ..."
# (restrict the pipeline), PORT=8000, IMAGE=sctcna:latest.

VENV      ?= .venv
PYTHON    ?= python3.12
PY        := $(VENV)/bin/python
APP_DATA  ?= data/app
PROCESSED ?= data/processed
RAW       ?= data/raw
DATASETS  ?=
PORT      ?= 8000
IMAGE     ?= sctcna:latest

# uv is the documented installer (see requirements.txt); fall back to venv + pip.
UV := $(shell command -v uv 2>/dev/null)
ifneq ($(UV),)
VENV_CREATE = $(UV) venv --python 3.12 $(VENV)
PIP_INSTALL = $(UV) pip install --python $(PY)
else
VENV_CREATE = $(PYTHON) -m venv $(VENV)
PIP_INSTALL = $(PY) -m pip install
endif

.PHONY: venv install data pipeline aggregate app test test-e2e lint format docker-build docker-run

venv:
	test -d $(VENV) || $(VENV_CREATE)

install: venv
	$(PIP_INSTALL) -r requirements.txt
	$(PY) -m playwright install chromium

data:
	$(PY) scripts/download_data.py --all --convert --raw $(RAW)

pipeline:
	$(PY) -m sctcna run --raw $(RAW) --out $(PROCESSED) $(if $(DATASETS),--datasets $(DATASETS),)

aggregate:
	$(PY) -m sctcna aggregate --raw $(RAW) --processed $(PROCESSED) --app $(APP_DATA)

app:
	SCTCNA_DATA_DIR=$(APP_DATA) $(VENV)/bin/uvicorn app.main:app --reload --host 127.0.0.1 --port $(PORT)

test:
	$(VENV)/bin/pytest -m "not e2e"

test-e2e:
	$(VENV)/bin/pytest -m e2e

lint:
	$(VENV)/bin/ruff check .
	$(VENV)/bin/ruff format --check .

format:
	$(VENV)/bin/ruff format .
	$(VENV)/bin/ruff check --fix .

docker-build:
	docker build -t $(IMAGE) .

docker-run:
	docker run --rm -p 7860:7860 $(IMAGE)
