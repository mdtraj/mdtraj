import mdtraj as md
from mdtraj.testing import get_fn
from sklearn.externals.joblib import load, dump

def test_1():
    # https://github.com/rmcgibbo/mdtraj/issues/438
    traj = md.load(get_fn('frame0.h5'))
    filenames = dump(traj, 'temp')
    print('saving %s' % ' '.join(filenames))
    traj = load('temp', mmap_mode='r')
    md.rmsd(traj, traj, 0) # BOOM

    for fn in filenames:
        os.unlink(fn)

if __name__ == '__main__':
    test_1()
