sudo apt-get update
sudo apt-get install -qq -y g++ gfortran valgrind csh
sudo apt-get install -qq -y g++-multilib gcc-multilib
wget http://repo.continuum.io/miniconda/Miniconda-3.6.0-Linux-x86_64.sh
bash Miniconda-3.6.0-Linux-x86_64.sh -b
PIP_ARGS="-U"

export PATH=$HOME/miniconda/bin:$PATH

conda update --yes conda
conda create --yes -n ${python} --file devtools/ci/requirements-conda-${python}.txt
source activate $python
$HOME/miniconda/envs/${python}/bin/pip install $PIP_ARGS -r devtools/ci/requirements-${python}.txt


# Install SPARTA+ for NMR chemical shift predicition
MDTRAJ_DIR=`pwd`
mkdir -p $HOME/external
cd $HOME/external
# wget http://spin.niddk.nih.gov/bax/software/SPARTA+/sparta+.tar.Z
wget https://s3-us-west-1.amazonaws.com/mdtraj-travis-ci-cached-files/sparta%2B.tar.Z
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

# -------------------

# Install shiftx2 for NMR chemical shift predicition
MDTRAJ_DIR=`pwd`
mkdir -p $HOME/external
cd $HOME/external
#wget http://www.shiftx2.ca/download/shiftx2-v107-linux-20120106.tgz
wget https://s3-us-west-1.amazonaws.com/mdtraj-travis-ci-cached-files/shiftx2-v107-linux-20120106.tgz
REFERENCE_MD5="4d3b23d77e773aa321af2a01ed04199a"
RECEIVED_MD5=$(md5sum shiftx2-v107-linux-20120106.tgz | cut -d " " -f 1)
if [ $REFERENCE_MD5 != $RECEIVED_MD5 ]; then
    echo "shiftx2-v107-linux-20120106.tgz md5 mismatch"
    exit 1
fi

tar -xzf shiftx2-v107-linux-20120106.tgz
cd shiftx2-v107-linux/
make
export PATH=`pwd`:$PATH
# go back to the original directory we were in
cd $MDTRAJ_DIR

# -------------------

# Install ppm for NMR chemical shift predicition
MDTRAJ_DIR=`pwd`
mkdir -p $HOME/external
cd $HOME/external
wget https://s3-us-west-1.amazonaws.com/mdtraj-travis-ci-cached-files/ppm_linux_64.exe
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
