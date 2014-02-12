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
