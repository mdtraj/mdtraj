sudo apt-get install -qq -y g++ gfortran valgrind
wget http://repo.continuum.io/miniconda/Miniconda-3.0.5-Linux-x86_64.sh
bash Miniconda-3.0.5-Linux-x86_64.sh -b

export PATH=$HOME/miniconda/bin:$PATH

conda create --yes -n ${python} python=${python} numpy scipy pytables cython \
    scipy pandas pip netcdf4 yaml nose

export PATH=$HOME/miniconda/envs/$python/bin:$PATH

PYTHON_VERSION=`python -c 'import sys; print("%d.%d" % sys.version_info[:2])'`
pip install $PIP_ARGS -r tools/ci/requirements-${PYTHON_VERSION}.txt