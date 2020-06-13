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
data = urlopen(URL + '/versions.json').read().decode()
versions = json.loads(data)

existing = [item['version'] for item in versions]
already_exists = version.short_version in existing

if not already_exists:
    # new release so all the others are now old
    for i in range(len(versions)):
        versions[i]['latest'] = False

    versions.append({
        'version': version.short_version,
        'display': version.short_version,
        'url': "{base}/{version}".format(base=URL, version=version.short_version),
        'latest': True})

with open("docs/_deploy/versions.json", 'w') as versionf:
    json.dump(versions, versionf, indent=2)

