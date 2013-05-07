"""
This file works at a frame-by-frame level.
"""
import numpy as np
import scipy.optimize

class Transformation():
    def __init__(self, rotation, translation):
        self.rotation = rotation
        self.translation = translation

    @classmethod
    def align(cls, xyz1, xyz2):
        """Returns a Transformation object to align xyz1 onto xyz2.

        Parameters
        ----------
        xyz : ndarray, shape = (n_atoms, 3)
            xyz coordinates of a SINGLE frame

        Returns
        -------
        R : ndarray
            rotation
        """    
        translation, rotation = find_translation_and_rotation(xyz1, xyz2)
        return cls(rotation, translation)    

    def transform(self, xyz):
        """Apply transformation to xyz."""
        mu1 = xyz.mean(0)
        return xyz.dot(self.rotation) + self.translation - mu1.dot(self.rotation)

def rmsd(xyz1, xyz2):
    """Returns a Transformation object to align xyz1 onto xyz2.

    Parameters
    ----------
    xyz : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a SINGLE frame

    Returns
    -------
    R : ndarray
        rotation
    """    
    T = Transformation.align(xyz1, xyz2)
    xyz1_prime = T.transform(xyz1)
    delta = xyz1_prime - xyz2
    rms = (delta ** 2.0).sum(1).mean() ** 0.5
    return rms

def find_translation_and_rotation(xyz1, xyz2):
    """Returns a Transformation object to align xyz1 onto xyz2.

    Parameters
    ----------
    xyz : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a SINGLE frame

    Returns
    -------
    R : ndarray
        rotation
    """    
    assert(xyz1.shape[1] == 3)
    assert(xyz1.shape == xyz2.shape)

    mu1 = xyz1.mean(0)
    mu2 = xyz2.mean(0)
    
    translation = mu2

    xyz1 = xyz1 - mu1
    xyz2 = xyz2 - mu2
    
    correlation_matrix = np.dot(np.transpose(xyz1), xyz2)
    V, S, W_tr = np.linalg.svd(correlation_matrix)
    is_reflection = (np.linalg.det(V) * np.linalg.det(W_tr)) < 0.0
    if is_reflection:
        V[:,-1] = -V[:,-1]
    rotation = np.dot(V, W_tr)
    
    return translation, rotation


"""
RMSD calculation, via quaterion-based characteristic polynomial in
pure python/numpy.

Reference
---------
Rapid calculation of RMSDs using a quaternion-based characteristic polynomial.
Acta Crystallogr A 61(4):478-480.


"""
import numpy as np

def _center(conformation):
    """Center and typecheck the conformation"""

    conformation = np.asarray(conformation)
    if not conformation.ndim == 2:
        raise ValueError('conformation must be two dimensional')
    _, three = conformation.shape
    if not three == 3:
        raise ValueError('conformation second dimension must be 3')

    centroid = np.mean(conformation, axis=0)
    centered = conformation - centroid
    return centered


def rmsd_qcp(conformation1, conformation2):
    """Compute the RMSD with Theobald's quaterion-based characteristic
    polynomial
    
    Rapid calculation of RMSDs using a quaternion-based characteristic polynomial.
    Acta Crystallogr A 61(4):478-480.
    
    Parameters
    ----------
    conformation1 : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the first conformation
    conformation2 : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the second conformation

    Returns
    -------
    rmsd : float
        The root-mean square deviation after alignment between the two pointsets
    """
    # center and typecheck the conformations
    A = _center(conformation1)
    B = _center(conformation2)
    if not A.shape[0] == B.shape[0]:
        raise ValueError('conformation1 and conformation2 must have same number of atoms')
    n_atoms = len(A)

    #the inner product of the structures A and B
    G_A = np.einsum('ij,ij', A, A)
    G_B = np.einsum('ij,ij', B, B)
    #print 'GA', G_A, np.trace(np.dot(A.T, A))
    #print 'GB', G_B, np.trace(np.dot(B.T, B))

    # M is the inner product of the matrices A and B
    M = np.dot(B.T, A)

    # unpack the elements
    Sxx, Sxy, Sxz = M[0, :]
    Syx, Syy, Syz = M[1, :]
    Szx, Szy, Szz = M[2, :]

    # do some intermediate computations to assemble the characteristic
    # polynomial
    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz

    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz

    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    # two of the coefficients
    C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C1 = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    # the other coefficient
    C0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 \
        + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) \
        + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz)) \
        + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz)) \
        + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz)) \
        + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz))

    E0 = (G_A + G_B) / 2.0
    f = lambda x: x ** 4.0 + C2 * x ** 2. + C1 * x + C0
    df = lambda x: 4 * x ** 3.0 + 2 * C2 * x + C1
    max_eigenvalue = scipy.optimize.newton(f, E0, df)        
    rmsd = np.sqrt(np.abs(2.0 * (E0 - max_eigenvalue) / n_atoms))
    return rmsd
