How to release
===================

Pre-release + Github
--------------------
- Update the `docs/whatsnew.rst` document. Use [this tool](https://github.com/rmcgibbo/gh-util),
  which should list all of the PRs that have been merged since the laster release.
- Update the version number in `setup.py`, change `ISRELEASED` to `True`
- Update the date in `docs/whatsnew.rst` and add a blurb about the release.
- Commit to master, and [tag](https://github.com/mdtraj/mdtraj/releases) the
  release on github.
- Make sure the docs were built and pushed to the correct place. 

PyPI
----
The next step is to add the release to the python package index.

- Run `git clean -fdx` to clean the source directory.
- Create the cannoncal "sdist" (source distribution) using `python setup.py sdist --formats=gztar,zip`.
- Inspect the sdist files (they're placed in `dist/`), and make sure they look right.
  You can try installing them into your environment with pip, unzipping or untaring them, etc.
- Once you're satisfied that the sdist is correct, push the source to PyPI using
  `twine upload [path to sdist files]`. This requires being registered on PyPI as a owner or maintainer
  of the project.

Immediately after creating the sdist
------------------------------------
- Update the version number in `setup.py` to `1.(y+1).0.dev0` per PEP440;
  change `ISRELEASED` to `False`.
- Add a new section in `docs/whatsnew.rst` and mark it "(Development)".
- Commit to master.

It's important that the version which is tagged on github for the release be
the (only) one with the ISRELEASED flag in setup.py set to true.

Point releases
--------------

If you want to release 1.y.1 after 1.y.0, create a branch named 1.y, apply
the fix, and follow the normal release procedure. If you need updated 
docs, change `travis.yml` to deploy on branch 1.y instead of master.

Conda
-----
- File a PR against [omnia-md/conda-recipes](https://github.com/omnia-md/conda-recipes)
- Update the recipe's version string and source URL to pull the new
  sdist from PyPI.
- Compare the dependency lists in meta.yaml to those in the travis tests.
  They should be the same!
- Travis and Appveyor will then build binary conda packages.
- Encourage people to test the "rc" release candidate packages with
  `conda install -c omnia/label/rc mdtraj`
- Move `rc` packages to `main` label via anaconda.org or posting an issue to
  omnia-md/conda-recipes.

Finally
-------
- Close the 1.y milestone. Make sure there is a 1.(y+1) milestone.
- Update this document with things you've learned!

Wheels
------
PyPI hosts *wheels*, pre-compiled binary packages, like conda packages, for OS X and
Windows. (At the time of this writing, they are still ironing out issues w.r.t.
linux.) To create and upload wheels, download the sdist and unpack the (or check out
the exact tag from git), and run `python setup.py bdist_wheel`.

For example, to build wheels for Python 2.7, 3.4 and 3.5 on OS X, I ran
```
conda env remove -y -n _build
versions=("2.7" "3.4" "3.5")
for v in "${versions[@]}"; do
    conda create -y -n _build python=$v numpy cython
    source activate _build
    python setup.py bdist_wheel
    source deactivate
    conda env remove -y -n _build
done
```
Then, if these all look good, you can upload them to PyPI with twine, as was done with the
sdist.


Docs Building & Hosting
=======================

After a travis build succeeds, the docs are built with sphinx and pushed to
the mdtraj.org amazon s3 account. The credentials for that account are stored,
encrypted, in the .travis.yml file.
(http://docs.travis-ci.com/user/build-configuration/#Secure-environment-variables)

As of version 0.7-dev (0.8 release) multiple versions of the docs are hosted
online. When a build happens on a version with ISRELEASED==False, it's put into
the "latest" folder on the S3 bucket. If ISRELEASED==True, it's put into a
subfolder with the name of the short release. The relevant logic is in
`tools/travis-ci/push-docs-to-s3.py`.

Tools License
=============
Copyright (c) 2012-2015 Stanford University and the Authors
All rights reserved.

Redistribution and use of all files in this folder (devtools) and (../.travis.yml,
../basesetup.py, ../setup.py) files in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
