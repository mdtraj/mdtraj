if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi

if [[ "$CONDA_PY" == "27" && "$CONDA_NPY" == "19" ]]; then
    # Create the docs and push them to S3

    conda create -yq -n py27 python=2.7
    source activate py27
    conda install -yq `conda build devtools/conda-recipe --output`
    pip install numpydoc s3cmd msmb_theme==0.3.1
    conda install -yq `cat docs/requirements.txt | xargs`

    (cd docs && make html)
    python devtools/travis-ci/push-docs-to-s3.py
    python devtools/travis-ci/update-versions.py
fi


