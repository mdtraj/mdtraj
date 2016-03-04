#!/bin/bash

# Install SPARTA+ for NMR chemical shift predicition
MDTRAJ_DIR=`pwd`
mkdir -p $HOME/external
cd $HOME/external
wget http://mdtraj.org/travis-ci-cache/sparta+.tar.Z
REFERENCE_MD5="d4293336254f5696221db0edcc57cfed"
RECEIVED_MD5=$(md5sum sparta+.tar.Z | cut -d " " -f 1)
if [ $REFERENCE_MD5 != $RECEIVED_MD5 ]; then
    echo "sparta+.tar.Z md5 mismatch"
    exit 1
fi

tar -xzf sparta+.tar.Z
cd SPARTA+
csh ./install.com
export SPARTAP_DIR=`pwd`
export SPARTA_DIR=`pwd`
export PATH=`pwd`/bin:$PATH
# go back to the original directory we were in
cd $MDTRAJ_DIR

