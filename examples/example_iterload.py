"""Out-of-core calculations with :func:`md.iterload`
"""
import numpy as np
np.set_printoptions(threshold=50)
import mdtraj as md

# Sometimes your molecular dynamics trajectory files are too large to
# fit into memory. This can make analysis a burden, because you have to
# be very aware of the size of various objects. This can be a challenge
# in python because of the language's automatic memory management.

# Fortunately, python provides the iterator protocol that can help
# us out here. We can "stream through" a trajectory, without loading
# the entire thing into memory at all. Instead, we'll process it in
# chunks.

# For the purpose of this example, we'll use a short trajectory
# that's included with MDTraj for testing purposes. When you use this
# recipe yourself, you probably will want to point your code to your
# own trajectory file

import mdtraj.testing
traj_filename = mdtraj.testing.get_fn('frame0.h5')

# First, if you only want a single frame of a trajectory, there's no reason
# to load up the whole thing. :func:`md.load_frame` can load up a single
# frame for you. Let's get the first one.
first_frame = md.load_frame(traj_filename, 0)
print first_frame

# Using :func:`md.iterload`, you can iterate through chunks
# of the trajectory. If you don't retain a reference to the chunk
# as you iterate through, then the python garbage collector can recycle
# the memory.

rmsds = []
for chunk in md.iterload(traj_filename, chunk=100):
    rmsds.append(md.rmsd(chunk, first_frame))
    print chunk
    print chunk.time

# Now, we've calculated all of the rmsds chunk by chunk, and we
# can take a look at them.

rmsds = np.concatenate(rmsds)
print rmsds
print np.max(rmsds), np.argmax(rmsds)
