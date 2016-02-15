
echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


if [[ "$CONDA_PY" == "27" -a "$CONDA_NPY" == "19" ]]; then
    # Create the docs and push them to S3
    # -----------------------------------

    conda create -n py27 python=2.7
    source activate py27
    conda install --yes `conda build devtools/conda-recipe --output`
    pip install numpydoc s3cmd msmb_theme==0.3.1
    conda install --yes `cat docs/requirements.txt | xargs`

    conda list -e

    (cd docs && make html)
    python devtools/travis-ci/push-docs-to-s3.py
    python devtools/travis-ci/update-versions.py
fi


