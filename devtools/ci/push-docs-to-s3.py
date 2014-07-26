import os
import boto
from boto.s3.key import Key
import mdtraj

# The secret key is available as a secure environment variable
# on travis-ci to push the build documentation to Amazon S3.
AWS_ACCESS_KEY_ID = os.environ['AWS_ACCESS_KEY_ID']
AWS_SECRET_ACCESS_KEY = os.environ['AWS_SECRET_ACCESS_KEY']
BUCKET_NAME = 'mdtraj.org'

bucket_name = AWS_ACCESS_KEY_ID.lower() + '-' + BUCKET_NAME
conn = boto.connect_s3(AWS_ACCESS_KEY_ID,
            AWS_SECRET_ACCESS_KEY)
bucket = conn.create_bucket(BUCKET_NAME)

root = 'docs/_build'
for dirpath, dirnames, filenames in os.walk(root):
    for filename in filenames:
        fn = os.path.join(dirpath, filename)
        print 'Uploading', fn, '...'
        k = Key(bucket)
        if not mdtraj.version.release:
            prefix = 'latest'
        else:
            prefix = mdtraj.version.short_version
        k.key =  os.path.join(prefix, os.path.relpath(fn, root))
        k.set_contents_from_filename(fn)
