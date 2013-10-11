# Compile wheels for dependencies in a vagrant VM. This reduces the time tests take on travis-ci

This only needs to be done when we want to change the version of dependencies. This is not
necessary very often

Instructions
------------
You need to get virtualbox and vagrant

To build all the wheels, we do it in a fresh precise64 virtual machine. This takes about
an hour. The wheels will be placed in the wheelhouse/ directory.
```
$ vagrant up
$ vagrant provision
```

Then you can push the wheelhouse to the web with `push-wheels-to-s3.py`

