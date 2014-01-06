"""
Registry for trajectory file formats, so that the appropriate file
object and loader can be resolved based on the filename extension.
"""

import os.path


class _FormatRegistry(object):
    """Registry for trajectory file objects.
    
    Example
    -------
    >>> @_FormatRegistry.register_loader('.xyz')
    >>> def load_xyz(filename):
        return Trajectory(...)
    >>>
    >>> print _FormatRegistry.loaders['.xyz']
    <function load_xyz at 0x1004a15f0>
    """
    loaders = {}
    fileobjects = {}

    @classmethod
    def register_loader(cls, extension):
        def decorator(f):
            cls.loaders[extension] = f
            return f
        return decorator

    @classmethod
    def register_fileobject(cls, extension):
        def decorator(f):
            cls.fileobjects[extension] = f
            return f
        return decorator

# Make a single instance of this class, and then
# get rid of the class object. This should be
# treated as a singleton
_FormatRegistry = _FormatRegistry()
