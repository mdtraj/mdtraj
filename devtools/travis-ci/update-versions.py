from __future__ import print_function
import os
import pip
import json
from tempfile import NamedTemporaryFile
import subprocess
from mdtraj import version
from six.moves.urllib.request import urlopen
if not any(d.project_name == 's3cmd' for d in pip.get_installed_distributions()):
    raise ImportError('The s3cmd pacakge is required. try $ pip install s3cmd')


URL = 'http://mdtraj.org/'
BUCKET_NAME = 'mdtraj.org'

if not version.release:
    print("This is not a release.")
    exit(0)


versions = json.load(urlopen(URL + '/versions.json'))

# new release so all the others are now old
for i in range(len(versions)):
    versions[i]['latest'] = False

versions.append({
    'version': version.short_version,
    'url': URL + '/' + str(version.short_version),
    'latest': True})

# The secret key is available as a secure environment variable
# on travis-ci to push the build documentation to Amazon S3.
with NamedTemporaryFile('w') as config, NamedTemporaryFile('w') as v:
    config.write('''[default]
access_key = {AWS_ACCESS_KEY_ID}
secret_key = {AWS_SECRET_ACCESS_KEY}
'''.format(**os.environ))
    json.dump(versions, v)
    config.flush()
    v.flush()

    template = ('s3cmd --config {config} '
                'put {vfile} s3://{bucket}/versions.json')
    cmd = template.format(
            config=config.name,
            vfile=v.name,
            bucket=BUCKET_NAME)
    subprocess.call(cmd.split())
