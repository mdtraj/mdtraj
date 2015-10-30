import os
import warnings
from IPython.display import display, Javascript
try:
    from notebook.nbextensions import install_nbextension
except ImportError:
    from IPython.html.nbextensions import install_nbextension

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import resource_filename

__all__ = ['enable_notebook']

_REQUIRE_CONFIG = Javascript('''
require.config({
    paths: {
        'three': '//cdnjs.cloudflare.com/ajax/libs/three.js/r68/three.min',
        'iview' : '/nbextensions/mdtraj/iview',
        'surface' : '/nbextensions/mdtraj/surface.min',
        'objexporter' : '/nbextensions/mdtraj/objexporter',
        'filesaver' : '/nbextensions/mdtraj/filesaver',
        'jqueryui': '//ajax.googleapis.com/ajax/libs/jqueryui/1.11.1/jquery-ui.min',
        'contextmenu': '/nbextensions/mdtraj/context',
    },
    shim: {
        three: {
            exports: 'THREE'
        },
        iview: {
            deps: ['three', 'surface'],
            exports: 'iview'
        },
        surface: {
            exports: 'ProteinSurface'
        },
        objexporter: {
            deps: ['three'],
            exports: 'THREE.OBJExporter'
        },
        jqueryui: {
            exports: "$"
        },
        contextmenu: {
            deps: ['jquery'],
            exports: "context"
        },
    },
});
''',
css  = ['//lab.jakiestfu.com/contextjs/context.standalone.css']
)

def enable_notebook():
    """Enable IPython notebook widgets to be displayed.

    This function should be called before using TrajectoryWidget.
    """
    display(_REQUIRE_CONFIG)
    staticdir =  resource_filename('mdtraj', os.path.join('html', 'static'))
    install_nbextension(staticdir, destination='mdtraj', user=True)
