#!/usr/bin/env bash
if [[ ${MDTRAJ_TEST_PEP517} = "YES" ]]; then
    echo "testing for PEP517, enable pip installing dependencies"
    export PIP_NO_INDEX=false
    export PIP_NO_DEPENDENCIES=false
    export PIP_IGNORE_INSTALLED=false
fi

env | sort
pip install .
