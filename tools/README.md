Developer Notes / Tools
=======================

Assorted notes for developers.

How to do a release
-------------------
- Update the whatsnew.rst document. Use the github view that shows all the commits to master since the last release to write it.
- Update the version number in `setup.py`, change `ISRELEASED` to `True`
- Commit to master, and [tag](https://github.com/rmcgibbo/mdtraj/releases) the release on github
- To push the source to PyPI, use `python setup.py sdist --formats=gztar,zip upload`
- Conda binaries need to built separately on each platform (`conda build mdtraj; binstar upload <path to .tar.bz2>`)
- Make an annoucement on github / email
- After tagging the release, make a NEW commit that changes `ISRELEASED` back to `False` in `setup.py`


It's important that the version which is tagged on github for the release be
the one with the ISRELEASED flag in setup.py set to true.


Docs Building & Hosting
-----------------------

After a travis build succeeds, the docs are built with sphinx and pushed to the mdtraj.org amazon s3 account. The credentials for that account are stored, encrypted, in the .travis.yml file. (http://docs.travis-ci.com/user/build-configuration/#Secure-environment-variables)

As of version 0.7-dev (0.8 release) multiple versions of the docs are hosted online. When a build happens on a version with ISRELEASED==False, it's put into the "latest" folder on the S3 bucket. If ISRELEASED==True, it's put into a subfolder with the name of the short release. The relevant logic is in `tools/ci/push-docs-to-s3.py`.

In order for the select bar at the bottom of the docs that toggles between versions to work, these folders MUST match up with the tag names on github. This is because the list of items to put in that little dropdown menu is dynamically determined from the github API in the browser. This is the only way I could think of to make sure the old docs have a link to the latest version. The logic that populates the version dropdown menu in the browser is in

`docs/themes/sphinx_rtd_theme-0.1.5/sphinx_rtd_theme/versions.html`

Specifically note that it goes to https://api.github.com/repos/rmcgibbo/mdtraj/releases, and uses the `tag_names` to build the links. So these much line up with the prefix of `mdtraj.version.short_version` used in `tools/ci/push-docs-to-s3.py` for the links not to break.