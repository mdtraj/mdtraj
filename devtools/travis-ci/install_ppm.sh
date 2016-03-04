#!/bin/bash

# Install ppm for NMR chemical shift predicition
MDTRAJ_DIR=`pwd`
mkdir -p $HOME/external
cd $HOME/external
wget http://mdtraj.org/travis-ci-cache/ppm_linux_64.exe
chmod a+x ppm_linux_64.exe
REFERENCE_MD5="f3cb5681bd2769cdcfc77fe17c563ee4"
RECEIVED_MD5=$(md5sum ppm_linux_64.exe | cut -d " " -f 1)
if [ $REFERENCE_MD5 != $RECEIVED_MD5 ]; then
    echo "ppm_linux_64.exe md5 mismatch"
    exit 1
fi

export PATH=`pwd`:$PATH
# go back to the original directory we were in
cd $MDTRAJ_DIR
