import xtc
import numpy as np
xyz = np.zeros((2, 22, 3), dtype=np.float32)

reader = xtc.XTCReader('/home/robert/msmbuilder/Tutorial/XTC/RUN00/frame0.xtc', 0, 0, 1)

reader.set_outarray(xyz)
reader.next()
reader.next()

print xyz
