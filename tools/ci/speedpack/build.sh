#!/bin/bash

# This script is meant to run on a mint precise64 VM.
# The generated wheel files should be compatible
# with travis-ci as of 07/2013.
#
# Runtime can be up to an hour or more.

echo "Building wheels..."

# print a trace for everything; RTFM
set -x

ls /refq
apt-get update
bash /reqf/apt_install.sh
pip install virtualenv
pip install -I git+https://github.com/pypa/pip@$pip_commit#egg=pip
pip install -I -U setuptools
pip install wheel

export PYTHONIOENCODING='utf-8'
export VIRTUALENV_DISTRIBUTE=0

function get_site_pkgs_dir() {
    python$1 -c 'import distutils; print(distutils.sysconfig.get_python_lib())'
}


function create_wheel() {
    local pip_args="$1"
    local wheelhouse="$2"
    local n="$3"
    local pyver="$4"

    local site_pkgs_dir="$(get_site_pkgs_dir $pyver)"


    if [[ "$n" == *statsmodels* ]]; then
        pip wheel $pip_args --wheel-dir=$wheelhouse $n && \
        pip install $pip_args --no-index $n && \
        rm -Rf $site_pkgs_dir
    else
        pip wheel $pip_args --wheel-dir=$wheelhouse $n
        pip install $pip_args --no-index $n
    fi
}


function generate_wheels() {
    # get the requirements file
    local reqfile="$1"

    # get the python version
    local TAG=$(echo $reqfile |  grep -Po "(\d\.?[\d\-](_\w+)?)")

    # base dir for wheel dirs
    local WHEELSTREET=/wheelhouse
    local WHEELHOUSE="$WHEELSTREET/$TAG"

    local PY_VER="${TAG:0:3}"
    local PY_MAJOR="${PY_VER:0:1}"
    local PIP_ARGS="--use-wheel --find-links=$WHEELHOUSE --download-cache /tmp"

    # install the python version if not installed
    apt-get install python$PY_VER python$PY_VER-dev -y

    # create a new virtualenv
    rm -Rf /tmp/venv
    virtualenv -p python$PY_VER /tmp/venv
    source /tmp/venv/bin/activate

    # install pip setuptools
    pip install -I --download-cache /tmp 'git+https://github.com/pypa/pip@42102e9d#egg=pip'
    pip install -I -U --download-cache /tmp setuptools
    pip install -I --download-cache /tmp wheel

    # make the dir if it doesn't exist
    mkdir -p $WHEELHOUSE

    # put the requirements file in the wheelhouse
    cp $reqfile $WHEELHOUSE

    # install and build the wheels
    cat $reqfile | while read N; do
        create_wheel "$PIP_ARGS" "$WHEELHOUSE" "$N" "$PY_VER"
    done
}


for reqfile in $(ls -1 /reqf/requirements-*.*); do
    generate_wheels "$reqfile"
done
