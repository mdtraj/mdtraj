import bz2
import gzip
import io
import os
from io import StringIO


def open_maybe_zipped(filename, mode, force_overwrite=True):
    """Open a file in text (not binary) mode, transparently handling
    .gz or .bz2 compresssion, with utf-8 encoding.

    Parameters
    ----------
    filename : str
        Path to file. Compression will be automatically detected if
        the filename ends in .gz or .bz2.
    mode : {'r', 'w'}
        Mode in which to open file
    force_overwrite : bool, default=True
        If 'w', should we overwrite the file if something with `filename`
        already exists?

    Returns
    -------
    handle : file
        Open file handle.
    """
    _, extension = os.path.splitext(str(filename).lower())
    if mode == "r":
        if extension == ".gz":
            with gzip.GzipFile(filename, "r") as gz_f:
                return StringIO(gz_f.read().decode("utf-8"))
        elif extension == ".bz2":
            with bz2.BZ2File(filename, "r") as bz2_f:
                return StringIO(bz2_f.read().decode("utf-8"))
        else:
            return open(filename)
    elif mode == "w":
        if os.path.exists(filename) and not force_overwrite:
            raise OSError('"%s" already exists' % filename)
        if extension == ".gz":
            binary_fh = gzip.GzipFile(filename, "wb")
            return io.TextIOWrapper(binary_fh, encoding="utf-8")
        elif extension == ".bz2":
            binary_fh = bz2.BZ2File(filename, "wb")
            return io.TextIOWrapper(binary_fh, encoding="utf-8")
        else:
            return open(filename, "w")
    else:
        raise ValueError('Invalid mode "%s"' % mode)
