import base64
import json
from IPython.display import display, Javascript
from IPython.html.widgets import DOMWidget
from IPython.utils.traitlets import Unicode, CUnicode, Bool, Bytes, CInt, Any, List, Dict

# Global. Inject JS into notebook
WIDGET_JS = Javascript(filename='static/widget_trajectory.js')
REQUIRE_CONFIG = Javascript('''
require.config({
    paths: {
        'three': '//cdnjs.cloudflare.com/ajax/libs/three.js/r68/three.min',
        'three/trackball' : 'http://mrdoob.github.io/three.js/examples/js/controls/TrackballControls',
        'rmol' : 'http://rawgit.com/rmcgibbo/mdtraj/notebook/MDTraj/html/libs/RMol',
    },
    shim: {
        'three': {
            exports: 'THREE'
        },
        'three/trackball': {
            deps: ['three'],
            exports: 'THREE'
        },
    },
});
console.log("setting shims");
''')

def enable_notebook():
    display(WIDGET_JS)
    display(REQUIRE_CONFIG)

class TrajectoryWidget(DOMWidget):
    disabled = Bool(False, help="Enable or disable user changes.", sync=True)
    _view_name = Unicode('TrajectoryView', sync=True)

    frame = CInt()
    trajectory = Any()
    _topology = Dict(sync=True)
    _xyz = List(sync=True)

    def __init__(self, trajectory, frame=0):
        super(TrajectoryWidget, self).__init__()
        self.trajectory = trajectory
        self.frame = frame

    def _frame_changed(self, name, old, new):
        self._xyz = self.trajectory.xyz[new].tolist()

    def _trajectory_changed(self, name, old, new):
        self._xyz = new.xyz[self.frame].tolist()
        self._topology = new.topology.to_dict()

    

