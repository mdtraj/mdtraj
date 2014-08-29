import base64
import time
import json
from IPython.display import display, Javascript, HTML
from IPython.html.widgets import DOMWidget, IntSliderWidget, ContainerWidget
from IPython.utils.traitlets import Unicode, CUnicode, Bool, Bytes, CInt, Any, List, Dict

# Global. Inject JS into notebook
WIDGET_JS = Javascript(filename='static/widget_trajectory.js')
REQUIRE_CONFIG = Javascript('''
require.config({
    paths: {
        'three': '//cdnjs.cloudflare.com/ajax/libs/three.js/r68/three.min',
        'three/trackball' : 'http://mrdoob.github.io/three.js/examples/js/controls/TrackballControls',
        // 'rmol' : 'http://rawgit.com/rmcgibbo/mdtraj/notebook/MDTraj/html/libs/RMol',
        'rmol' : '//localhost:8000/libs/RMol',
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
    display(HTML('''
    <script src="//cdnjs.cloudflare.com/ajax/libs/three.js/r68/three.min.js"></script>
    <script src="//localhost:8000/libs/TrackballControls.js"></script>
    <script src="//localhost:8000/libs/RMol.js"></script>
    <script src="//localhost:8000/static/widget_trajectory.js"></script>
    <script src="//localhost:8000/static/widget_imagebutton.js"></script>
    '''))
    # display(WIDGET_JS)
    # display(REQUIRE_CONFIG)

class TrajectoryFrameWidget(DOMWidget):
    disabled = Bool(False, help="Enable or disable user changes.", sync=True)
    _view_name = Unicode('TrajectoryView', sync=True)

    frame = CInt()
    trajectory = Any()
    _topology = Dict(sync=True)
    _xyz = List(sync=True)

    def __init__(self, trajectory, frame=0):
        super(TrajectoryFrameWidget, self).__init__()
        self.trajectory = trajectory
        self.frame = frame

    def _frame_changed(self, name, old, new):
        self._xyz = self.trajectory.xyz[new].tolist()

    def _trajectory_changed(self, name, old, new):
        self._xyz = new.xyz[self.frame].tolist()
        self._topology = new.topology.to_dict()


def TrajectoryWidget(trajectory):
    viewer = TrajectoryFrameWidget(trajectory)
    slider = IntSliderWidget(min=0, max=trajectory.n_frames)

    def on_value_change(name, value):
        viewer.frame = value
    slider.on_trait_change(on_value_change, 'value')
    return ContainerWidget(children=(viewer, slider))

