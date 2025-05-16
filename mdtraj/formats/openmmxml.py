##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import numpy as np

from mdtraj.formats.registry import FormatRegistry

__all__ = ["load_xml"]


@FormatRegistry.register_loader(".xml")
def load_xml(filename, top=None):
    """Load a single conformation from an OpenMM XML file.

    The OpenMM serialized state XML format contains additional information that
    is not read by this method, including forces, energies, and velocities.
    Here, we just read the positions and the box vectors.

    Parameters
    ----------
    filename : path-like
        The path on disk to the XML file
    top : {str, Trajectory, Topology}
        The XML format does not contain topology information. Pass in either the
        path to a pdb file, a trajectory, or a topology to supply this information.

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.
    """
    import xml.etree.ElementTree as etree

    from mdtraj.core.trajectory import Trajectory, _parse_topology

    topology = _parse_topology(top)

    tree = etree.parse(filename)

    # get all of the positions from the XML into a list of tuples
    # then convert to a numpy array
    positions = []
    for position in tree.getroot().find("Positions"):
        positions.append(
            (
                float(position.attrib["x"]),
                float(position.attrib["y"]),
                float(position.attrib["z"]),
            ),
        )

    box = []
    vectors = tree.getroot().find("PeriodicBoxVectors")
    for name in ["A", "B", "C"]:
        box.append(
            (
                float(vectors.find(name).attrib["x"]),
                float(vectors.find(name).attrib["y"]),
                float(vectors.find(name).attrib["z"]),
            ),
        )

    traj = Trajectory(xyz=np.array(positions), topology=topology)
    traj.unitcell_vectors = np.array(box).reshape(1, 3, 3)

    return traj
