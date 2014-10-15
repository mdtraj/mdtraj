
Atom Selection
==============

Basics
------

In this example, we'll go over the basics of atom and reside selection
in MDTraj. First let's load up an example trajectory.

.. code:: python

    import mdtraj as md
    traj = md.load('ala2.h5')
    print traj

.. parsed-literal::

    <mdtraj.Trajectory with 100 frames, 22 atoms, 3 residues, without unitcells>


We can also more directly find out how many atoms or residues there are
by using ``traj.n_atoms`` and ``traj.n_residues``.

.. code:: python

    print 'How many atoms? %s' % traj.n_atoms
    print 'How many residues? %s' % traj.n_residues

.. parsed-literal::

    How many atoms? 22
    How many residues? 3


We can also manipulate the atom positions by working with ``traj.xyz``,
which is a NumPy array contain the xyz coordinated of each atom with
dimensions (n\_frames, n\_atoms, 3). Let's find the 3D coordinates of
the tenth atom in frame 5.

.. code:: python

    frame_idx = 4 # zero indexed frame number
    atom_idx = 9 # zero indexed atom index
    print 'Where is the fifth atom at the tenth frame?'
    print 'x: %s\ty: %s\tz: %s' % tuple(traj.xyz[frame_idx, atom_idx,:])

.. parsed-literal::

    Where is the fifth atom at the tenth frame?
    x: 0.697151	y: 0.92419	z: 0.872604


Topology Object
---------------

As mentioned previously in the introduction, every ``Trajectory`` object
contains a ``Topology``. The ``Topology`` of a ``Trajectory`` contains
all the connectivity information of your system and specific chain,
residue, and atom information.

.. code:: python

    topology = traj.topology
    print topology

.. parsed-literal::

    <mdtraj.Topology with 1 chains, 3 residues, 22 atoms, 21 bonds>


With the topology object we can select a certain ``atom``, or loop
through them all. (Note: everything is zero-indexed)

.. code:: python

    print 'Fifth atom: %s' % topology.atom(4)
    print 'All atoms: %s' % [atom for atom in topology.atoms]

.. parsed-literal::

    Fifth atom: ACE1-C
    All atoms: [ACE1-H1, ACE1-CH3, ACE1-H2, ACE1-H3, ACE1-C, ACE1-O, ALA2-N, ALA2-H, ALA2-CA, ALA2-HA, ALA2-CB, ALA2-HB1, ALA2-HB2, ALA2-HB3, ALA2-C, ALA2-O, NME3-N, NME3-H, NME3-C, NME3-H1, NME3-H2, NME3-H3]


The same goes for residues.

.. code:: python

    print 'Second residue: %s' % traj.topology.residue(1)
    print 'All residues: %s' % [residue for residue in traj.topology.residues]

.. parsed-literal::

    Second residue: ALA2
    All residues: [ACE1, ALA2, NME3]


Additionally, every ``atom`` and ``residue`` is also an object, and has
it's own set of properties. Here is a simple example that showcases just
a few.

.. code:: python

    atom = topology.atom(10)
    print '''Hi! I am the %sth atom, and my name is %s. 
    I am a %s atom with %s bonds. 
    I am part of an %s residue.''' % ( atom.index, atom.name, atom.element.name, atom.n_bonds, atom.residue.name)                                                                                                

.. parsed-literal::

    Hi! I am the 10th atom, and my name is CB. 
    I am a carbon atom with 4 bonds. 
    I am part of an ALA residue.


There are also more complex properties, like ``atom.is_sidechain`` or
``residue.is_protein``, which allow for more powerful selections.

Putting Everything Together
---------------------------

Hopefully, you can see how these properties can be combined with
Python's filtered list functionality. Let's say we want the indices of
all carbon atoms in the sidechains of our molecule. We could try
something like this.

.. code:: python

    print [atom.index for atom in topology.atoms if atom.element.symbol is 'C' and atom.is_sidechain]

.. parsed-literal::

    [1, 4, 8, 10, 14, 18]


Or maybe we want all even-indexed residues in the first chain (Although
this example only has the one chain....).

.. code:: python

    print [residue for residue in topology.chain(0).residues if mod(residue.index,2) == 0]

.. parsed-literal::

    [ACE1, NME3]


Atom Selection Language
-----------------------

If you're hesistant about programming filtered lists like the ones
above, MDTraj also features a rich atom selection language, similar to
that of PyMol and VMD. You can access it by using ``topology.select``.
Let's find all atoms in the last two residues.

.. code:: python

    print topology.select('resid 1 to 2')

.. parsed-literal::

    [ 6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21]


You can also do more complex operations. Here, we're looking for all
nitrogen atoms in the backbone.

.. code:: python

    print topology.select('name N and backbone')

.. parsed-literal::

    [ 6 16]


If you ever want to see the code that generates these results you can
use ``select_expression``, which will yield a string represention of the
atom selection code.

.. code:: python

    selection = topology.select_expression('name CA and resid 1 to 2')
    print selection

.. parsed-literal::

    [atom.index for atom in topology.atoms if ((atom.name == 'CA') and (1 <= atom.residue.index <= 2))]

