from __future__ import print_function

import os
import datetime

import boto
from boto.s3.key import Key
from jinja2 import Environment


# The secret key is available as a secure environment variable
# on travis-ci to push the build documentation to Amazon S3.
AWS_ACCESS_KEY_ID = os.environ['AWS_ACCESS_KEY_ID']
AWS_SECRET_ACCESS_KEY = os.environ['AWS_SECRET_ACCESS_KEY']
BUCKET_NAME = 'mdtraj-deps-wheelhouse'
ROOT = 'wheelhouse'

def main():
    bucket_name = AWS_ACCESS_KEY_ID.lower() + '-' + BUCKET_NAME
    conn = boto.connect_s3(AWS_ACCESS_KEY_ID,
                           AWS_SECRET_ACCESS_KEY)
    bucket = conn.create_bucket(BUCKET_NAME)

    for dirpath, dirnames, filenames in os.walk(ROOT):
        for filename in filenames:
            fn = os.path.join(dirpath, filename)
            print('Uploading', fn, '...')
            k = Key(bucket)
            k.key = os.path.relpath(fn, ROOT)
            k.set_contents_from_filename(fn)
        
        # Put up an index page
        k = Key(bucket)
        k.key = os.path.relpath(os.path.join(dirpath, 'index.html'), ROOT)
        k.metadata={'Content-Type': 'text/html'}
        k.set_contents_from_string(index_page(dirpath, filenames))


def index_page(dirpath, filenames):
    # Standard index page, copied from apache

    html = '''
    <html>
     <head>
      <title>Index of {{ dirpath }} </title>
     </head>
     <body>
    <h1>Index of {{ dirpath }} </h1>
    <table><tr><th><img src="/icons/blank.gif" alt="[ICO]"></th><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=M;O=A">Last modified</a></th><th><a href="?C=S;O=A">Size</a></th><th><a href="?C=D;O=A">Description</a></th></tr><tr><th colspan="5"><hr></th></tr>

    {% for fn in filenames %}
    <tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="{{ fn }}">{{ fn }} </a></td><td align="right"> {{ ts }} </td><td align="    right"> </td><td>&nbsp;</td></tr>
    {% endfor %}

    <tr><th colspan="5"><hr></th></tr>
    </table>
    </body></html>'''
    return Environment().from_string(html).render(dirpath=dirpath, filenames=filenames, ts=str(datetime.datetime.now()))


if __name__ == '__main__':
    main()
