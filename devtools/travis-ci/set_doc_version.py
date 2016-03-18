import os
import shutil
from mdtraj import version

if version.release:
    docversion = version.short_version
else:
    docversion = 'development'

os.mkdir("docs/_deploy")
shutil.copytree("docs/_build/html", "docs/_deploy/{docversion}"
                .format(docversion=docversion))
