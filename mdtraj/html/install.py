import os
import warnings
from IPython.display import display, Javascript
from IPython.html.nbextensions import install_nbextension
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import resource_filename

__all__ = ['enable_notebook']

_REQUIRE_CONFIG = Javascript('''
require.config({
    paths: {
        'three': '//cdnjs.cloudflare.com/ajax/libs/three.js/r68/three.min',
        'iview' : '/nbextensions/iview',
        'surface' : '/nbextensions/surface.min',
        'exporter' : '/nbextensions/objexporter',
        'filesaver' : '/nbextensions/filesaver',
        'jqueryui': '//ajax.googleapis.com/ajax/libs/jqueryui/1.11.1/jquery-ui.min',
        'contextmenu': '/nbextensions/context',
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
        exporter: {
            deps: ['three'],
            exports: 'THREE.OBJExporter'
        },
        jqueryui: {
            exports: "$"
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
    libs = ['iview.js','surface.min.js','objexporter.js','filesaver.js','context.js']
    fns = [resource_filename('mdtraj', os.path.join('html', 'static', f)) for f in libs]
    for fn in fns:
        install_nbextension(fn, verbose=0, user=True)
    display(_REQUIRE_CONFIG)

    widgets = ['widget_trajectory.js', 'widget_imagebutton.js']
    for fn in widgets:
        fn = resource_filename('mdtraj', os.path.join('html', 'static', fn))
        display(Javascript(filename=fn))
