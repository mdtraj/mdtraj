import numpy as np
import warnings

class TypeCastPerformanceWarning(RuntimeWarning):
    pass

def ensure_type(val, dtype, ndim, name, length=None, can_be_none=False, shape=None,
    warn_on_cast=True):
    """Ensure dtype and shape of an ndarray
    
    Parameters
    ----------
    val : {np.ndaraay, None}
        The array to check
    dtype : {nd.dtype, str}
        The dtype you'd like the array to have
    ndim : int
        The number of dimensions you'd like the array to have
    name : str
        name of the array. This is used when throwing exceptions.
    length : int, optional
        How long should the array be?
    can_be_none : bool
        Is `val=None` acceptable?
    shape : tuple, optional
        What should be shape of the array be
    warn_on_cast : bool, default=True
        Raise a warning when the dtypes don't match and a cast is done.
    
    Notes
    -----
    The returned value will always be C-contiguous.
        
    Returns
    -------
    typechecked_val : np.ndarray, None
        If `val=None` and `can_be_none=True`, then this will return None.
        Otherwise, it will return val (or a copy of val). If the dtype wasn't right,
        it'll be casted to the right shape. If the array was not C-contiguous, it'll
        be copied as well.
    
    """
    if can_be_none and val is None:
        return None

    if not isinstance(val, np.ndarray):
        raise TypeError(("%s must be numpy array. "
            " You supplied type %s" % (name, type(val))))
    
    if warn_on_cast and val.dtype != dtype:
        warnings.warn("Casting dtype=%s to %s" % (val.dtype, dtype),
            TypeCastPerformanceWarning)
        
    val = np.ascontiguousarray(val, dtype=dtype)
    if not val.ndim == ndim:
        raise ValueError(("%s must be ndim %s. "
            "You supplied %s" % (name, ndim, val.ndim)))

    if length is not None and len(val) != length:
        raise ValueError(("%s must be length %s. "
            "You supplied %s" % (name, length, len(val))))

    if shape is not None and val.shape != shape:
        raise ValueError(("%s must be shape %s. "
            "You supplied %s" % (name, shape, val.shape)))

    return val