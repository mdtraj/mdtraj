#!/usr/bin/env bash
if [[ ${MDTRAJ_TEST_PEP517} = "YES" ]]; then
    echo "testing for PEP517, enable pip installing dependencies"
    # since we want build time dependencies to be installed, we have to unset conda builds safety guards
    # hindering us from doing so. Then we build a mdtraj wheel, which we install with the old settings.
    _PIP_NO_INDEX=$PIP_NO_INDEX
    _PIP_NO_DEPENDENCIES=$PIP_NO_DEPENDENCIES
    _PIP_IGNORE_INSTALLED=$PIP_IGNORE_INSTALLED
    _PIP_NO_BUILD_ISOLATION=$PIP_NO_BUILD_ISOLATION

    unset PIP_NO_INDEX
    unset PIP_NO_DEPENDENCIES
    unset PIP_IGNORE_INSTALLED
    unset PIP_NO_BUILD_ISOLATION
    pip wheel .

    export PIP_NO_INDEX=$_PIP_NO_INDEX
    export PIP_NO_DEPENDENCIES=$_PIP_NO_DEPENDENCIES
    export PIP_IGNORE_INSTALLED=$_PIP_IGNORE_INSTALLED
    export PIP_NO_BUILD_ISOLATION=$_PIP_NO_BUILD_ISOLATION
    pip install mdtraj*.whl -v
else
    pip install .
fi
