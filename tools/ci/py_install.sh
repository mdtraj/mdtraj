
pip install -I git+https://github.com/pypa/pip@$pip_commit#egg=pip
pip install -I -U setuptools
pip install wheel

base_url=http://stanford.edu/~rmcgibbo/wheelhouse
PIP_ARGS=" -I --use-wheel --find-links=$base_url/${TRAVIS_PYTHON_VERSION}/"
pip install $PIP_ARGS -r requirements-${TRAVIS_PYTHON_VERSION}.txt

