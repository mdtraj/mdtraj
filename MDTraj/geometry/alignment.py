"""
This file works at a frame-by-frame level.
"""
import numpy as np

class Transformation():
    def __init__(self, rotation, translation):
        self.rotation = rotation
        self.translation = translation

    @classmethod
    def align(cls, xyz1, xyz2):
        """Returns a Transformation object to align xyz1 onto xyz2.

        Parameters
        ----------
        xyz : ndarray
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
    xyz : ndarray
        xyz coordinates of a SINGLE frame

    Returns
    -------
    R : ndarray
        rotation
    """    
    T = Transformation.align(xyz1, xyz2)
    xyz1_prime = T.transform(xyz1)
    delta = xyz1_prime - xyz2
    rms = (delta ** 2.0).mean() ** 0.5
    return rms

def find_translation_and_rotation(xyz1, xyz2):
    """Returns a Transformation object to align xyz1 onto xyz2.

    Parameters
    ----------
    xyz : ndarray
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
