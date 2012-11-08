from ctypes import *
import sys

dcdplugin=cdll.LoadLibrary("../dcdplugin_s.so")

# Wrapper for the timestep structure. 
class MolfileTimestep(Structure):
    _fields_ = [
        ("coords",POINTER(c_float)),
        ("velocities",POINTER(c_float)),
        ("A",c_float),
        ("B",c_float),
        ("C",c_float),
        ("alpha",c_float),
        ("beta",c_float),
        ("gamma",c_float),
        ("physical_time",c_double)
        ]

dcdplugin.open_dcd_read.argtypes=[c_char_p,c_char_p,POINTER(c_int)]
dcdplugin.open_dcd_read.restype=c_void_p

# 0 is OK,  -1 is EOF, else error
dcdplugin.read_next_timestep.argtypes=[c_void_p,c_int,POINTER(MolfileTimestep)]
dcdplugin.read_next_timestep.restype=c_int

dcdplugin.close_file_read.argtypes=[c_void_p]
dcdplugin.close_file_read.restype=None

# The following can also be exported
#  open_dcd_write
#  write_timestep
#  close_file_write

fn=sys.argv[1]
na=c_int(-1)


init_err=dcdplugin.vmdplugin_init()
assert init_err==0

v=dcdplugin.open_dcd_read(fn,"dcd",byref(na))
print "Number of atoms (66): %d"% na.value

Coords=c_float * (na.value*3)
xvec=Coords()

ts=MolfileTimestep()
ts.coords=xvec

frame_no=0

while True:
    ts_err=dcdplugin.read_next_timestep(v,na,byref(ts))
    assert ts_err>=-1
    # frame 0, index 0, X
    if (frame_no==0):
        print "Frame 0, index 0, X: %f. Delta: %f" %  (ts.coords[0],ts.coords[0]+0.73600000143051147)
    # frame 1, index 1, Z
    elif (frame_no==1):
        print "Frame 1, index 1, Z: %f. Delta: %f" %  (ts.coords[0],ts.coords[5]-2.872999906539917)
    if ts_err==-1:
        break
    frame_no=frame_no+1

print "Read %d frames" % frame_no
print "Physical time at end: %f" % ts.physical_time
print "Box at end: %f %f %f   %f %f %f" % (ts.A, ts.B, ts.C, ts.alpha, ts.beta, ts.gamma)



