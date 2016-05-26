# coding: utf-8
def call_read(n_frames):
    data = fh._read(n_frames, None)
    print (data.shape, data)
from mdtraj.testing import get_fn
from mdtraj.formats.tng import TNGTrajectoryFile
f=get_fn('tng_example.tng')
fh = TNGTrajectoryFile(f)
call_read(4)
call_read(4)
call_read(4)
call_read(4)

fh.close()
