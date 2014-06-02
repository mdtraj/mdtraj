coveralls

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" == "true" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "2.7 3.3 3.4" =~ "$python" ]]; then
    conda install binstar
    binstar -t $BINSTAR_TOKEN  upload -u omnia -p mdtraj-dev $HOME/miniconda/conda-bld/linux-64/mdtraj-dev-*
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi

# Create the docs and push them to S3
# -----------------------------------

# Install stuff for running the example IPython notebooks
sudo apt-get install -qq pandoc         # notebook -> rst
conda install --yes matplotlib scikit-learn sphinx boto ipython-notebook jinja2
pip install runipy==0.0.4                      # example notebooks

# Install OpenMM for a couple of the the examples
conda config --add channels http://conda.binstar.org/omnia
conda install --yes openmm

cd docs && make html && cd -
python tools/ci/push-docs-to-s3.py
