from __future__ import absolute_import
import base64
from IPython.display import display, Javascript
from IPython.html.widgets import DOMWidget, IntSliderWidget, ContainerWidget
from IPython.utils.traitlets import Unicode, Bool, Bytes, CInt, Any, List, Dict, Enum


class TrajectoryWidget(DOMWidget):
    disabled = Bool(False, help="Enable or disable user changes.", sync=True)
    _view_name = Unicode('TrajectoryView', sync=True)

    frame = CInt()
    trajectory = Any()
    topology = Dict(sync=True)
    coordinates = List(sync=True)
    
    # Display options
    height = CInt('300', help='Height in pixels', sync=True)
    width = CInt('300', help='Width in pixles', sync=True)
    color = Enum(['chainbow', 'ss', 'chain', 'polarity'], 'chainbow', sync=True)
    mainChain = Enum(['ribbon', 'thickRibbon', 'strand', 'chain',
                      'cylinderHelix', 'tube', 'bonds'],
                     'thickRibbon', sync=True)
    sideChains = Enum(['line', None], None, sync=True)

    def __init__(self, trajectory, frame=0, **kwargs):        
        super(TrajectoryWidget, self).__init__(**kwargs)
        self.trajectory = trajectory
        self.frame = frame

    def _frame_changed(self, name, old, new):
        self.coordinates = self.trajectory.xyz[new].tolist()

    def _trajectory_changed(self, name, old, new):
        self.coordinates = new.xyz[self.frame].tolist()
        self.topology = new.topology.to_dict()
