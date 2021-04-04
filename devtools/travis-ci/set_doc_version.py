import os
import shutil
from mdtraj import version

if version.release:
    docversion = version.short_version
else:
    docversion = 'development'

os.makedirs("docs/_deploy", exist_ok=True)
shutil.copytree("docs/_build/html", f"docs/_deploy/{docversion}")
