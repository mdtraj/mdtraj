import os
import pip
import tempfile
import subprocess
import mdtraj.version


BUCKET_NAME = 'mdtraj.org'
if not mdtraj.version.release:
    PREFIX = 'latest'
else:
    PREFIX = mdtraj.version.short_version

if not any(d.project_name == 's3cmd' for d in pip.get_installed_distributions()):
    raise ImportError('The s3cmd pacakge is required. try $ pip install s3cmd')

# The secret key is available as a secure environment variable
# on travis-ci to push the build documentation to Amazon S3.
with tempfile.NamedTemporaryFile('w') as f:
    f.write('''[default]
access_key = {AWS_ACCESS_KEY_ID}
secret_key = {AWS_SECRET_ACCESS_KEY}
'''.format(**os.environ))
    f.flush()

    template = ('s3cmd --config {config} '
                'sync doc/_build/ s3://{bucket}/{prefix}/')
    cmd = template.format(
            config=f.name,
            bucket=BUCKET_NAME,
            prefix=PREFIX)
    subprocess.call(cmd.split())
