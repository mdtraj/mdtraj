import json

try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
from mdtraj import version

if not version.release:
    print("This is not a release.")
    exit(0)

URL = 'http://www.mdtraj.org'
versions = json.load(urlopen(URL + '/versions.json'))

# new release so all the others are now old
for i in range(len(versions)):
    versions[i]['latest'] = False

versions.append({
    'version': version.short_version,
    'display': version.short_version,
    'url': "{base}/{version}".format(base=URL, version=version.short_version),
    'latest': True})

with open("doc/_deploy/versions.json", 'w') as versionf:
    json.dump(versions, versionf)

