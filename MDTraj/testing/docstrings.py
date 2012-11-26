from inspect import (isclass, ismodule, isfunction, ismethod, \
    getmembers, getdoc, getmodule, getargs, isbuiltin)
import docscrape
import sys

__all__ = ['DocStringFormatTester']

def DocStringFormatTester(module, error_on_none=False):
    """
    Short Descr

    Parameters
    ----------
    module : module
        Descr
    error_on_none : bool
        Throw an error if no docstring is defined

    """
    import types

    # these are the functions whose docstrings we make methods to check
    def accept(f):
        # the builtin one is specifically for cython functions, which get
        # called builtins even though they're not
        return isfunction(f) or ismethod(f) or isbuiltin(f)
    functions = filter(accept, walk(module))

    def format(f):
        if ismethod(f):
            return '.'.join([getmodule(f).__name__, f.im_class.__name__, f.__name__])
        if isfunction(f) or isbuiltin(f):
            return '.'.join([getmodule(f).__name__, f.__name__])
        if isclass(f):
            return f.__name__
        return 'Error'

    def check_docstring(self, f):
        doc = getdoc(f)
        if doc is None:
            if error_on_none:
                raise ValueError('no docstring for %s' % format(f))
        else:
            parsed = docscrape.NumpyDocString(doc)
            param_names = set([e[0] for e in parsed['Parameters']])
            
            if isbuiltin(f):
                # You can't get the arglist from a builtin, which
                # is how cython functions turn up
                
                # but you can, hackily, get the number of arguments it wants
                # by parseing the error hen you supply too many
                import re
                try:
                    f(*range(1000))
                except TypeError as e:
                    m = re.search('takes at most (\d+) positional arguments', e.message)
                    if not m:
                        return
                    n_args = int(m.group(1))
                
                if len(param_names) != n_args:
                    raise ValueError("In %s, number of arguments, %d, doesn't "
                        " match the length of the Parameters in the "
                        "docstring, %d" % (format(f), n_args, len(param_names)))
                return
                

            args = set(getargs(f.func_code).args)
            if ismethod(f):
                args.remove('self')

            if args != param_names:
                raise ValueError("In %s, arguments %s don't "
                    "match Parameters list %s" % (format(f),
                        list(args), list(param_names)))


    funcdict = {}
    for i, f in enumerate(functions):
        name = 'test_%s' % i
        func = types.FunctionType(check_docstring.func_code, globals(), name,
            (f,), check_docstring.func_closure)

        func.func_doc = 'NumpyDoc: ' + format(f)
        funcdict[name] = func

    return type('TestDoc', (), funcdict)


def ispackage(obj):
    """
    Short Descr

    Parameters
    ----------
    obj : module
        Descr
    """
    if ismodule(obj):
        return obj.__file__.endswith("__init__.pyc") or \
            obj.__file__.endswith("__init__.py")
    return False


def walk(module):
    """
    Short Descr

    Parameters
    ----------
    module : module
        Descr
    """
    assert ismodule(module)
    if ispackage(module):
        raise ValueError('No packages')

    def is_valid(obj):
        if getmodule(obj) == module:
            # cython specific stuff
            if module.__file__.endswith('.so'):
                if isbuiltin(obj):
                    return not obj.__name__.startswith('_')

            if ismethod(obj) or isfunction(obj):
                return not obj.__name__.startswith('_')
            if isclass(obj):
                return True
        return False


    instack = [v for k, v in getmembers(module) if is_valid(v)]
    outstack = []

    while True:
        try:
            item = instack.pop()
        except IndexError:
            break
        outstack.append(item)

        if isclass(item):
            instack.extend([v for k, v in getmembers(item) if is_valid(v)])

    return outstack


TestModule = DocStringFormatTester(sys.modules[__name__])
