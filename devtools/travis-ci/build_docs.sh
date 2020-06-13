#!/bin/bash

# Print each line, exit on error
set -ev


# Install doc requirements
conda config --set safety_checks disabled
conda config --add channels omnia
conda config --add channels conda-forge
conda install -yq --file docs/requirements.txt
pip install -U awscli msmb_theme==1.2.0
which gcc
which g++
python setup.py install

# Make docs
cd docs && make html && cd -

# Move the docs into a versioned subdirectory
python devtools/travis-ci/set_doc_version.py

# Prepare versions.json
python devtools/travis-ci/update_versions_json.py
