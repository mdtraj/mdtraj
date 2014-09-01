from __future__ import absolute_import
import base64
from itertools import groupby

import mdtraj as md

from IPython.display import display, Javascript
from IPython.html.widgets import DOMWidget, IntSliderWidget, ContainerWidget
from IPython.utils.traitlets import Unicode, Bool, Bytes, CInt, Any, List, Dict, Enum


class TrajectoryWidget(DOMWidget):
    """IPython notebook widget for displaying trajectories in the browser with WebGL

    Example
    -------
    # if the final line occurs at the end of an IPython notebook cell, the
    # resulting interactive widget will be displayed
    >>> t = md.load('trajectory.pdb')
    >>> from mdtraj.html import enable_notebook, TrajectoryWidget
    >>> enable_notebook()
    >>> widget = TrajectoryWidget(t)
    >>> widget

    Attributes
    ----------
    frame : int
        Index of the frame to display.
    height : int
        Height, in pixels, of the display window
    width : int
        Width, in pixels, of the display window
    color : {'chainbow', 'ss', 'chain', 'polarity'}
        Color scheme used for the protein display
    mainChain : {'ribbon', 'thickRibbon', 'strand', 'chain', 'cylinderHelix', 'tube', 'bonds'}
        Drawing scheme for the main protein chain
    sideChains : {'line', None}
        Drawing scheme for the sidechains

    Notes
    -----
    All of the attributes listed above are synced with the browser's widget.
    Modifying these attributes, after the widget is constructed, will cause
    the widget to update *live*. They can also be set at widget construction
    time as keyword arguments to ``__init__``.

    See Also
    --------
    enable_notebook() : Executing this function before using the widget is
        required to load the required browser-side libraries
    """
    disabled = Bool(False, help="Enable or disable user changes.", sync=True)

    # Name of the javascript class which this widget syncs against on the
    # browser side. To work correctly, this javascript class has to be
    # registered and loaded in the browser before this widget is constructed
    # (that's what enable_notebook() does)
    _view_name = Unicode('TrajectoryView', sync=True)

    frame = CInt(0, help='Which frame from the trajectory to display')
    trajectory = Any()

    # The essence of the IPython interactive widget API on the python side is
    # that by declaring traitlets with sync=True, these variables are
    # automatically synced by the IPython runtime between this class in Python
    # and the browser-side model. Changes to these attributes are propagated
    # automatically to the browser (and changes on the browser side can trigger
    # events on this class too, although we're not using that feature).
    _topology = Dict(sync=True)
    _secondaryStructure = Dict(sync=True)
    _coordinates = List(sync=True)

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
        self._secondaryStructure = self._computeSecondaryStructure()

    def _trajectory_changed(self, name, old, new):
        self._coordinates = new.xyz[self.frame].tolist()
        self._topology = new.topology.to_dict()
        self._secondaryStructure = self._computeSecondaryStructure()

    def _computeSecondaryStructure(self):
        SS_MAP = {'C': 'coil', 'H': 'helix', 'E': 'sheet'}

        top = self.trajectory.topology
        dssp = md.compute_dssp(self.trajectory[self.frame])[0]
        result = {}

        # iterate over the (rindx, ss) pairs in enumerate(dssp),
        # and use itertools to group them into streaks by contiguous
        # chain and ss.
        keyfunc = lambda ir : (top.residue(ir[0]).chain, ir[1])
        for (chain, ss), grouper in groupby(enumerate(dssp), keyfunc):
            # rindxs is a list of residue indices in this contiguous run
            rindxs = [g[0] for g in grouper]
            for r in rindxs:
                # add entry for each atom in the residue
                for a in top.residue(r).atoms:
                    result[a.index] = {
                        'ss': SS_MAP[ss],
                        'ssbegin': (r==rindxs[0] and ss in {'H', 'E'}),
                        'ssend': (r==rindxs[-1] and ss in {'H', 'E'})}
        return result
