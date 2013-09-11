import os
import boto
from boto.s3.key import Key
import mdtraj

# Fill these in - you get them when you sign up for S3
AWS_ACCESS_KEY_ID = 'AKIAJN6O6JPBACO27BLQ'
AWS_SECRET_ACCESS_KEY = os.environ['AWS_SECRET_ACCESS_KEY']

bucket_name = AWS_ACCESS_KEY_ID.lower() + '-mdtraj'
conn = boto.connect_s3(AWS_ACCESS_KEY_ID,
            AWS_SECRET_ACCESS_KEY)
bucket = conn.create_bucket('mdtraj')

for dirpath, dirnames, filenames in os.walk('../mdtraj-docs'):
    for filename in filenames:
        fn = os.path.join(dirpath, filename)
        print 'Uploading', fn, '...'
        k = Key(bucket)
        k.key = os.path.relpath(fn, '../mdtraj-docs')
        k.set_contents_from_filename(fn)
