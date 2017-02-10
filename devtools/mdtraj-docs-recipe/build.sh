#!/bin/bash

# TODO: Make conda packages
pip install sphinx_rtd_theme==0.1.9 msmb_theme==1.2.0
sphinx-build -b html -d $SRC_DIR/doctrees $SRC_DIR/docs $PREFIX/share/mdtraj-docs
