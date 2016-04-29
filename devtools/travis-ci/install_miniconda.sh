#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]];   then MINICONDA=Miniconda3-latest-MacOSX-x86_64.sh; fi
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then MINICONDA=Miniconda3-latest-Linux-x86_64.sh;  fi

MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b
rm -f $MINICONDA

export PATH=$HOME/miniconda3/bin:$PATH

conda update -yq conda
conda install -yq conda-build jinja2
