import os
import shutil
import sys

from mdtraj import version

if version.release:
    docversion = version.short_version
else:
    docversion = 'development'

os.mkdir("docs/_deploy")
shutil.copytree("{prefix}/share/mdtraj-docs".format(prefix=sys.prefix),
                "docs/_deploy/{docversion}".format(docversion=docversion))
