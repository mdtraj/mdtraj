#!/usr/bin/env bash
if [[ -n MDTRAJ_TEST_PEP517 ]]; then
    echo "testing for PEP517, enable pip installing dependencies"
    export PIP_NO_INDEX=false
    export PIP_NO_DEPENDENCIES=false
    #export PIP_IGNORE_INSTALLED=false
fi

env | sort
pip install -vv --no-deps .
