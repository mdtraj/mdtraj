import contextlib
import os
import shutil
import tempfile
import time

__all__ = ["timing", "enter_temp_directory"]


class timing:
    """A timing context manager

    Examples
    --------
    >>> long_function = lambda : None
    >>> with timing('long_function'):
    ...     long_function()
    long_function: 0.000 seconds
    """

    def __init__(self, name="block"):
        self.name = name
        self.time = 0
        self.start = None
        self.end = None

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, ty, val, tb):
        self.end = time.time()
        self.time = self.end - self.start
        print(f"{self.name}: {self.time:0.3f} seconds")
        return False


@contextlib.contextmanager
def enter_temp_directory():
    """Create and enter a temporary directory; used as context manager."""
    temp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(temp_dir)
    yield
    os.chdir(cwd)
    shutil.rmtree(temp_dir)
