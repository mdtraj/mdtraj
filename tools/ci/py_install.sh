set -x
pip_commit=42102e9deaea99db08b681d06906c2945f6f95e2
PYTHON_VERSION=`python -c 'import sys; print("%d.%d" % sys.version_info[:2])'`
pip install -I git+https://github.com/pypa/pip@$pip_commit#egg=pip
pip install -I -U setuptools
pip install wheel


base_url=http://stanford.edu/~rmcgibbo/wheelhouse
PIP_ARGS=" -I --use-wheel --find-links=$base_url/${PYTHON_VERSION}/"
pip install $PIP_ARGS -r tools/ci/requirements-${PYTHON_VERSION}.txt

