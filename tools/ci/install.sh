sudo apt-get install -qq -y g++ gfortran valgrind
wget http://repo.continuum.io/miniconda/Miniconda-3.0.5-Linux-x86_64.sh
bash Miniconda-3.0.5-Linux-x86_64.sh -b

export PATH=$HOME/miniconda/bin:$PATH

conda create --yes -n ${python} python=${python} numpy scipy pytables cython \
    scipy pandas pip netcdf4 yaml nose

export PATH=$HOME/miniconda/envs/$python/bin:$PATH

PYTHON_VERSION=`python -c 'import sys; print("%d.%d" % sys.version_info[:2])'`
pip install $PIP_ARGS -r tools/ci/requirements-${PYTHON_VERSION}.txt


# Install SPARTA+ for NMR
mkdir $HOME/external
cd $HOME/external
wget http://spin.niddk.nih.gov/bax/software/SPARTA+/sparta+.tar.Z
REFERENCE_MD5="12a04ec45d9bd9e7974b218fe2353765"
RECEIVED_MD5=$(md5sum sparta+.tar.Z | cut -d " " -f 1)
if [$MD5_1 != $MD5_2] then
    echo "sparta+.tar.Z md5 mismatch"
    exit 1
fi

tar -xzvf sparta+.tar.Z
cd SPARTA+
./install.com
export SPARTAP_DIR=$HOME/external/SPARTA+
export SPARTA_DIR=$HOME/external/SPARTA+
export PATH=$HOME/external/SPARTA+:$PATH
