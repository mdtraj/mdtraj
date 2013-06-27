.. _building-docs:

Building the documentation
--------------------------

To build the documentation and push it to github, it's a little tricky.
Unfortunately, I'm not sure if we can use readthedocs, since they don't support
compiling c code which is necessary to get the docstrings for our cython modules.
First, you'll need to install sphinx and numpydoc ::

    pip install sphinx numpydoc
  
Next, make a new directory for the documentation that's not inside the mdtraj
repository. It should be "adjacent to it". From the mdtraj project root, run ::

    cd /path/to/mdtraj
    cd ..
    git clone git@github.com:rmcgibbo/mdtraj.git mdtraj-docs
    cd mdtraj-docs
    git checkout gh-pages
    
Now, go back to the docs subdirectory in the main repository. This will require
push access to the mdtraj repository. ::

    cd /path/to/mdtraj/docs
    make html
    make github
