coveralls

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


anaconda -t $ANACONDA_TOKEN  upload --force -u omnia -p mdtraj-dev $HOME/miniconda/conda-bld/linux-64/mdtraj-dev-*


if [[ "$python" == "2.7" ]]; then
    # Create the docs and push them to S3
    # -----------------------------------

    conda install --yes `conda build devtools/conda-recipe --output`
    pip install numpydoc s3cmd msmb_theme
    conda install --yes `cat docs/requirements.txt | xargs`

    conda list -e

    (cd docs && make html)
    python devtools/travis-ci/push-docs-to-s3.py
    python devtools/travis-ci/update-versions.py
fi


