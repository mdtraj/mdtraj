pip install -I git+https://github.com/pypa/pip@$pip_commit#egg=pip
pip install -I -U setuptools
pip install wheel
PYTHON_VERSION=`python -c 'import sys; print("%d.%d" % sys.version_info[:2])'`

base_url=http://stanford.edu/~rmcgibbo/wheelhouse
PIP_ARGS=" -I --use-wheel --find-links=$base_url/${PYTHON_VERSION}/"
pip install $PIP_ARGS -r requirements-${PYTHON_VERSION}.txt

