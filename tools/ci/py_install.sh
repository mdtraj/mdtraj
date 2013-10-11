set -x

$python -V

# Remove any old distribute or setuptools
sudo rm -rf `$python -c 'import distutils.sysconfig; print(distutils.sysconfig.get_python_lib())'`/setuptools*
sudo rm -rf `$python -c 'import distutils.sysconfig; print(distutils.sysconfig.get_python_lib())'`/distribute*

# Get  a new setuptools
wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py -O - | sudo $python
# Get a bootstrap pip
sudo easy_install pip

# Install Dependencies
# as of pip 1.4rc2, wheel files are still being broken regularly, this is a
# known good commit. should revert to pypi when a final release is out
pip_commit=42102e9deaea99db08b681d06906c2945f6f95e2
PYTHON_VERSION=`$python -c 'import sys; print("%d.%d" % sys.version_info[:2])'`
sudo pip install -I git+https://github.com/pypa/pip@$pip_commit#egg=pip
sudo pip install -I -U wheel


base_url=http://stanford.edu/~rmcgibbo/wheelhouse
PIP_ARGS=" -I --use-wheel --find-links=$base_url/${PYTHON_VERSION}/"
sudo pip install $PIP_ARGS -r tools/ci/requirements-${PYTHON_VERSION}.txt

