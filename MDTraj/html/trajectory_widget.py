from __future__ import absolute_import
import base64
from itertools import groupby

import mdtraj as md

from IPython.display import display, Javascript
from IPython.html.widgets import DOMWidget, IntSliderWidget, ContainerWidget
from IPython.utils.traitlets import (Unicode, Bool, Bytes, CInt, Any,
                                     Dict, Enum)

__all__ = ['TrajectoryView']


class TrajectoryView(DOMWidget):
    """IPython notebook widget for displaying trajectories in the browser with WebGL

    Example
    -------
    # if the final line occurs at the end of an IPython notebook cell, the
    # resulting interactive widget will be displayed
    >>> t = md.load('trajectory.pdb')
    >>> from mdtraj.html import enable_notebook, TrajectoryView
    >>> enable_notebook()
    >>> widget = TrajectoryView(t)
    >>> widget

    Attributes
    ----------
    camera : {'perspective', 'orthographic'}
        Camera mode (default='perspective')
    background : {'black', 'grey', 'white'}
        Background color (default='black')
    colorBy : {'spectrum', 'chain', 'secondary structure', 'residue',
               'polarity', 'atom'}
        Color scheme (default='white')
    primaryStructure : {'lines', 'stick', 'ball & stick','sphere', 'nothing'}
        Drawing method for the primary structure (default='nothing')
    secondaryStructure = {'ribbon', 'strand', 'cylinder & plate', 'C alpha trace', 'nothing'}
        Drawing method for secondary structure. (default='cylinder & plate')
    surfaceRepresentation = {'Van der Waals surface', 'solvent excluded surface',
                                'solvent accessible surface', 'molecular surface', 'nothing'}
        Drawing method for surface representation. (default='nothing')

    Notes
    -----
    All of the attributes listed above are synced with the browser's widget.
    Modifying these attributes, after the widget is constructed, will cause
    the widget to update *live*. They can also be set at widget construction
    time as keyword arguments to ``__init__``.

    The viewer WebGL viewer used, iview, is documented in [1].

    References
    ----------
    ..[1] Li, Hongjian, et al. "iview: an interactive WebGL visualizer for
          protein-ligand complex." BMC Bioinformatics 15.1 (2014): 56.

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
    _frameData = Dict(sync=True)

    # Display options
    camera = Enum(['perspective', 'orthographic'], 'perspective', sync=True)
    background = Enum(['black', 'grey', 'white'], 'white', sync=True)
    colorBy = Enum(['spectrum', 'chain', 'secondary structure', 'residue',
                    'polarity', 'atom'], 'spectrum', sync=True)
    primaryStructure = Enum(['lines', 'stick', 'ball & stick', 'sphere',
                             'nothing'], 'nothing', sync=True)
    secondaryStructure = Enum(['ribbon', 'strand', 'cylinder & plate',
                               'C alpha trace', 'nothing'], 'cylinder & plate',
                               sync=True)
    surfaceRepresentation = Enum(['Van der Waals surface','solvent excluded surface', 
                    'solvent accessible surface', 'molecular surface',
                    'nothing'], 'nothing', sync=True)
    
    def __init__(self, trajectory, frame=0, **kwargs):
        super(TrajectoryView, self).__init__(**kwargs)
        self.trajectory = trajectory
        self.frame = frame

    def _frame_changed(self, name, old, new):
        """Automatically called by the traitlet system when self.frame is modified"""
        self._update_frame_data()

    def _trajectory_changed(self, name, old, new):
        """Automatically called by the traitlet system when self.trajectory is modified"""
        self._topology = self._computeTopology()
        self._update_frame_data()


    def _update_frame_data(self):
        self._frameData = {
            'coordinates' : self.trajectory.xyz[self.frame].tolist(),
            'secondaryStructure' : self._computeSecondaryStructure()
        }

    def _computeSecondaryStructure(self):
        """Compute the secondary structure of the selected frame and
        format it for the browser
        """
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
                        'ssbegin': (r==rindxs[0] and ss in set(['H', 'E'])),
                        'ssend': (r==rindxs[-1] and ss in set(['H', 'E']))}
        return result

    def _computeTopology(self):
        """Extract the topology and format it for the browser. iview has a
        particular format for storing topology-based information, and
        for simplicity and hack-ability the best place to do the
        conversion is here in python.
        """
        # TODO(rmcgibbo). Document this schema. It needs to match with what's
        # going on inside iview.loadTopology on the browser side.

        
        atoms = {}

        # these should be mutually exclusive. you're only in one of
        # these categories
        peptideIndices = []
        waterIndices = []
        ionIndices = []
        ligandIndices = []

        bondIndices = []
        calphaIndices = []

        hetIndices = []

        for atom in self.trajectory.topology.atoms:
            atoms[atom.index] = {
                'alt': ' ',
                'b' : 0,
                'chain' : atom.residue.chain.index,
                'elem': atom.element.symbol if atom.element is not None else 'X',
                'insc' : ' ',
                'name' : atom.name,
                'resi' : atom.residue.index,
                'resn' : atom.residue.name,
                'serial' : atom.index,
                'ss' : None,
                'coord' : None,
                'bonds' : [],
            }
            if atom.name == 'CA':
                calphaIndices.append(atom.index)

            if atom.residue.is_water:
                waterIndices.append(atom.index)
            elif not atom.residue.is_protein:
                ligandIndices.append(atom.index)
            else:
                peptideIndices.append(atom.index)

        for ai, aj in self.trajectory.topology.bonds:
            bondIndices.append((ai.index, aj.index))

        return {'atoms': atoms,
                'bondIndices' : bondIndices,
                'ionIndices': ionIndices,
                'calphaIndices': calphaIndices,
                'hetIndices': hetIndices,
                'peptideIndices': peptideIndices,
                'ligandIndices' : ligandIndices,
                'waterIndices' : waterIndices}
