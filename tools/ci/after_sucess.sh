PYTHON_VERSION=`python -c 'import sys; print("%d.%d" % sys.version_info[:2])'`

if [[ "$TRAVIS_PULL_REQUEST" == "true" ]]; then 
    echo "This is a pull request. No deployment will be done."; exit 0
fi

if [[ "$PYTHON_VERSION" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${PYTHON_VERSION}"; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


sudo apt-get install -qq python-sphinx       # for building documentation
sudo apt-get install -qq python-boto         # for interacting with S3
sudo pip install ipython
cd docs && make html && cd -
python tools/ci/push-docs-to-s3.py
