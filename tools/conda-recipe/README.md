This is a recipe for building the current development package into a conda
binary.

The installation on travis-ci is done by building the conda package, installing
it, running the tests, and then if successful pushing the package to binstar
(and the docs to AWS S3). The binstar auth token is an encrpyed environment
variable generated using:

binstar auth -n mdtraj-travis -o omnia --max-age 22896000 -c --scopes api:write

and then saved in the environment variable BINSTAR_TOKEN.


