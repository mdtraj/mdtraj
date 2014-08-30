from __future__ import absolute_import
import base64
from IPython.display import display, Javascript
from IPython.html.widgets import DOMWidget, IntSliderWidget, ContainerWidget
from IPython.utils.traitlets import Unicode, Bool, Bytes, CInt, Any, List, Dict


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

