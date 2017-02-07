#!/bin/bash

# Print each line, exit on error
set -ev

# Install the built package
conda create --yes -n docenv python=$CONDA_PY
source activate docenv
conda install -yq --use-local mdtraj mdtraj-docs

# Move the docs into a versioned subdirectory
python devtools/travis-ci/set_doc_version.py

# Prepare versions.json
python devtools/travis-ci/update_versions_json.py
