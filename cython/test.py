import xtc

filename = '/Users/rmcgibbo/local/msmbuilder/Tutorial/XTC/RUN00/frame0.xtc'

# o = xtc.XTCReader(filename, chunk_len=50)
# 
# print o.load()

print len(xtc.load_xtc(filename, chunk_len=50)[3])

# print o.box
# 
# import msmbuilder.xtc
# 
# config = msmbuilder.xtc.XTCReader(filename).next()
# print config.box
