Developer Notes / Tools / License
=================================

Assorted notes for developers.

How to do a release
-------------------
- Update the `docs/whatsnew.rst` document. Use the github view that shows all the
  commits to master since the last release to write it.
- Update the version number in `setup.py`, change `ISRELEASED` to `True`
- Update the version number in `devtools/conda-recipe/meta.yaml`
- Commit to master, and [tag](https://github.com/mdtraj/mdtraj/releases) the
  release on github
- To push the source to PyPI, use `python setup.py sdist --formats=gztar,zip upload`.
  **Before pushing to PyPI, run `git clean -x -f -d` to make sure that nothing
  else gets added to the source release by accident**
- Update `meta.yaml` in the mdtraj/mdtraj-conda-binaries repository. Commit
  to master. Travis and Appveyor will build binary conda packages.
- After tagging the release, make a NEW commit that changes `ISRELEASED` back
  to `False` in `setup.py`


It's important that the version which is tagged on github for the release be
the one with the ISRELEASED flag in setup.py set to true.


Docs Building & Hosting
-----------------------

After a travis build succeeds, the docs are built with sphinx and pushed to
the mdtraj.org amazon s3 account. The credentials for that account are stored,
encrypted, in the .travis.yml file.
(http://docs.travis-ci.com/user/build-configuration/#Secure-environment-variables)

As of version 0.7-dev (0.8 release) multiple versions of the docs are hosted
online. When a build happens on a version with ISRELEASED==False, it's put into
the "latest" folder on the S3 bucket. If ISRELEASED==True, it's put into a
subfolder with the name of the short release. The relevant logic is in
`tools/travis-ci/push-docs-to-s3.py`.

In order for the select bar at the bottom of the docs that toggles between
versions to work, these folders MUST match up with the tag names on github.
This is because the list of items to put in that little dropdown menu is
dynamically determined from the github API in the browser. This is the only
way I could think of to make sure the old docs have a link to the latest
version. The logic that populates the version dropdown menu in the browser is in

`docs/themes/sphinx_rtd_theme-0.1.5/sphinx_rtd_theme/versions.html`

Specifically note that it goes to https://api.github.com/repos/mdtraj/mdtraj/releases,
and uses the `tag_names` to build the links. So these must line up with the
prefix of `mdtraj.version.short_version` used in `tools/travis-ci/push-docs-to-s3.py`
for the links not to break.

License
-------
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
