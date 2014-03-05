.. _HDF5FormatSpec:

MDTraj HDF5 Format Specification
================================

Overview
--------

This is the specification for a molecular dynamics trajectory file
format based on HDF5 which is supported by the MDTraj package.

Why HDF5?
---------

Design Goals
~~~~~~~~~~~~

In storing MD trajectory data for the for purposes including very large
scale analysis, there are a few design goals. (1) The trajectories should
be small and space efficient on disk. (2) The trajectories should be fast to
write and fast to read. (3) The data format should support flexible read
options. For instance, random access to different frames in the
trajectory should be possible. It should be possible to easily query the
dimensions of the trajectory (``n_frames``, ``n_atoms``, etc) without
loading the file into memory. It should be possible to load only every
n-th frame, or to directly only a subset of the atoms with limited
memory overhead. (5) The trajectory format should be easily extensible
in a backward compatible manner. For instance, it should be possible to
add new arrays/fields like the potential energy or the topology without
breaking backwards compatibility.

Other Formats
~~~~~~~~~~~~~

Currently, MDTraj is able to read and write trajectories in DCD, XTC,
TRR, BINPOS, and AMBER NetCDF formats, in addition to HDF5. This
presents an opportunity to compare these formats and see how they fit
our design goals. The most space efficient is XTC, because it uses 16 bit
fixed precision encoding. For some reason, the XTC read times are quite
slow though. DCD is fast to read and write, but relatively inflexible.
NetCDF is fast and flexible. BINPOS and MDCRD are garbage -- they're
neither fast, small, nor flexible.

What's lacking?
~~~~~~~~~~~~~~~

Of the formats we currently have, AMBER NetCDF is the best, in that it
it satisfies all of the design goals except for the first. But the
trajectories are twice as big on disk as XTC, which is really quite
unfortunate. For dealing with large data sets, size matters. So let's
define a HDF5 standard that has the benefits of AMBER NetCDF and the
benefits of XTC mixed together. We'll use an extensible data format
(HDF5), we'll provide options for lossy and lossless compression, and
we'll store the topology inside the trajectory, so that a single
trajectory file always contains the information needed to understand
(and visualize) the system.

Details
-------

This specification is heavily influenced by the AMBER NetCDF
`standard <http://ambermd.org/netcdf/nctraj.html>`__. Significant
portions of the text are copied verbatim.

Encoding
~~~~~~~~

-  Files will be encoded in HDF5, `a data model, library, and file
   format for storing and managing
   data <http://www.hdfgroup.org/HDF5/>`__ produced at NCSA.
-  Arrays may be encoded with zlib compression.
-  Libraries implementing this standard may, at their desecration, round
   the data to an appropriate number of significant digits, which can
   significantly enhance zlib compression ratios.

Recommended file extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  The ``.h5`` extension will be preferred, although not required.

Global Attributes
~~~~~~~~~~~~~~~~~

-  **Conventions (required)** Contents of this attribute are a comma or
   space delimited list of tokens representing all of the conventions to
   which the file conforms. Creators shall include the string Pande as
   one of the tokens in this list. In the usual case, where the file
   conforms only to this convention, the value of the attribute will
   simply be \`\`Pande''. Readers may fail if this attribute is not
   present or none of the tokens in the list are Pande. Optionally, if
   the reader does not expect HDF5 files other than those conforming to
   the Pande convention, it may emit a warning and attempt to read the
   file even when the Conventions attribute is missing.
-  **ConventionVersion (required)** Contents are a string representation
   of the version number of this convention. Future revisions of this
   convention having the same version number may include definitions of
   additional variables, dimensions or attributes, but are guaranteed to
   have no incompatible changes to variables, dimensions or attributes
   specified in previous revisions. Creators shall set this attribute to
   "1.1". If this attribute is present and has a value other than "1.1",
   readers may fail or may emit a warning and continue. It is expected
   that the version of this convention will change rarely, if ever.
-  **application (optional)** If the creator is part of a suite of
   programs or modules, this attribute shall be set to the name of the
   suite.
-  **program (required)** Creators shall set this attribute to the name
   of the creating program or module.
-  **programVersion (required)** Creators shall set this attribute to
   the preferred textual formatting of the current version number of the
   creating program or module.
-  **title (optional)** Creators may set use this attribute to represent
   a user-defined title for the data represented in the file. Absence of
   a title may be indicated by omitting the attribute or by including it
   with an empty string value.
-  **randomState (optional), ASCII-encoded string** Creators may
   optionally describe the state of their internal random number
   generators at the start of their simulation. The semantics of this
   string are specific to the MD code and are not specified by this
   standard.
-  **forcefield (optional), ASCII-encoded string** For data from a
   molecular dynamics simulation, creators may optionally describe the
   Hamiltonian used. This should be a short, human readable string, like
   "AMBER99sbildn".
-  **reference (optional), ASCII-encoded string** Creators may
   optionally specify published a reference that documents the program
   or parameters used to generate the data. The reference should be
   listed in a simple, human readable format. Multiple references may be
   listed simply by separating the references with a human readable
   delimiter within the string, like a newline.

Arrays
~~~~~~

-  **coordinates (required) shape=(n\_frames, n\_atoms, 3),
   type=Float32, units="nanometers"**. This variable shall contain the
   Cartesian coordinates of the specified particle for the specified.
-  **time (optional), shape=(n\_frames), dtype=Float32,
   units="picoseconds"** When coordinates on the frame dimension have a
   temporal sequence (e.g. they form a molecular dynamics trajectory),
   creators shall define this dimension and write a float for each frame
   coordinate representing the simulated time value in picoseconds
   associated with the frame. Time zero is arbitrary, but typically will
   correspond to the start of the simulation. When the file stores a
   collection of conformations having no temporal sequence, creators
   shall omit this variable.
-  **cell\_lengths (optional), shape=(n\_frames, 3, 3), dtype=Float32,
   units="nanometers"** When the data in the coordinates variable come
   from a simulation with periodic boundaries, creators shall include
   this variable. his variable shall represent the lengths (a,b,c) of
   the unit cell for each frame. The edge with length a lies along the x
   axis; the edge with length b lies in the x-y plane. The origin (point
   of invariance under scaling) of the unit cell is defined as (0,0,0).
   If the simulation has one or two dimensional periodicity, then the
   length(s) corresponding to spatial dimensions in which there is no
   periodicity shall be set to zero.
-  **cell\_angles shape=(n\_frames, 3, 3), dtype=Float32,
   units="degrees"** Creators shall include this variable if and only if
   they include the cell\_lengths variable. This variable shall
   represent the angles (, , ) defining the unit cell for each frame.
   defines the angle between the b and c vectors, defines the angle
   between the a and c vectors and defines the angle between the a and b
   vectors. Angles that are undefined due to less than three dimensional
   periodicity shall be set to zero.
-  **velocities (optional), shape=(n\_frames, n\_atoms, 3),
   type=Float32, units="nanometers/picosecond"** When the velocities
   variable is present, it shall represent the cartesian components of
   the velocity for the specified particle and frame. It is recognized
   that due to the nature of commonly used integrators in molecular
   dynamics, it may not be possible for the creator to write a set of
   velocities corresponding to exactly the same point in time as defined
   by the time variable and represented in the coordinates variable. In
   such cases, the creator shall write a set of velocities from the
   nearest point in time to that represented by the specified frame.
-  **kineticEnergy (optional), shape=(n\_frames), type=Float32,
   units="kJ/mol"** Creators may optionally specify the kinetic energy
   of the system at each frame.
-  **potentialEnergy (optional), shape=(n\_frames), type=Float32,
   units="kJ/mol"** Creators may optionally specify the potential energy
   of the system at each frame.
-  **temperature (optional), shape=(n\_frames), type=Float32,
   units="Kelvin"** Creators may optionally specify the temperature of
   the system at each frame.
-  **lambda (optional), shape=(n\_frames), type=Floa32 units=""** For
   describing an alchemical free energy simulation, a creator may
   optionally notate each frame in the simulation with a value of
   lambda.
-  **constraints (optional), shape=(n\_constraints, 3),
   type=CompoundType(int, int, float) units=[None, None, "nanometers"]**
   Creators may optionally describe any constraints applied to the bond
   lengths. ``constraints`` shall be a compound-type table (referred to
   a table as opposed to an array in the pytables documentation), such
   that the first two entries are the indices of the two atoms involved
   in the constant, and the final entry is the distance those atoms are
   constrained to.
-  **topology (optional, but highly recommended), shape=(1,
   length\_as\_needed) type=string** For protein systems, creators shall
   describe the topology of the system in ASCII encoded JSON. The format
   for the topology definition is described in the topology subsection
   of this document. The JSON string encoding the topology shall be
   stored as the sole row in an array of strings.


Array Metadata
~~~~~~~~~~~~~~

-  For arrays that contain naturally unitted numbers (which is all of
   them except for 'topology'), creators shall explicitly declare their
   units. The unit system of length=nanometers, time=picoseconds,
   mass=daltons, temperature=Kelvin, energy=kJ/mol, force=kJ/mol/nm
   shall be used everywhere. For angles, degrees shall be used. The
   units shall be set as an "attribute", on the array, under the key
   "units", within the parlance of HDF5. It shall be a string.

-  For arrays that contain numbers which have been rounded to a certain
   number of significant digits, creators shall declare the number of
   significant digits by setting the "least\_significant\_digit"
   attribue, which should be a positive integer.


Extended Arrays
~~~~~~~~~~~~~~~
Creators may extend this format by adding new arrays. Arrays containing
per-atom and per-frame data that naturally possesses physical units should
declare those units explicitly in the array attributes. Readers should be
flexible, ignoring the presence of arrays that they are not equiped to handle.


Topology
--------

Rational
~~~~~~~~

It is our experience that not having the topology stored in the same
file as the the trajectory's coordinate data is a pain -- it's just really
inconvenient. And generally, the trajectories are long enough that it
doesn't take up much incremental storage space to store the topology in
there too. The topology is not that complicated.

Format
~~~~~~

The topology will be stored in JSON. The JSON will then be serialized as
a string and stored in the HDF5 file with an ASCII encoding.

The topology stores a hierarchical description of the chains, residues,
and atoms in the system. Each chain is associated with an ``index`` and
a list of residues. Each residue is associated with a ``name``, an
``index``, a ``resSeq`` index (not zero-indexed), and a list of ``atom``\ s. 
Each ``atom`` is associated with a
``name``, an ``element``, and an ``index``. All of the indicies should
be zero-based.

The ``name`` of a residue is not strictly proscribed, but should
generally follow PDB 3.0 nomenclature. The ``element`` of an atom
shall be one of the one or two letter element abbreviations from the
periodic table. The ``name`` of an atom shall indicate some information
about the type of the atom beyond just its element, such as 'CA' for
the alpha carbom, 'HG' for a gamma hydrogen, etc. This format
does not specify exactly what atom names are allowed -- creators should
follow the conventions from the forcefield they are using.

In addition to the chains, the topology shall also contain a list of the
bonds. The bonds shall be a list of length-2 lists of integers, where
the integers refer to the index of the two ``atoms`` that are bonded.

Example
~~~~~~~

The following shows the topology of alanine dipeptide in this format.
Since it's JSON, the whitespace is optional and just for readability.

::

    {'bonds': [[4, 1],
               [4, 5],
               [1, 0],
               [1, 2],
               [1, 3],
               [4, 6],
               [14, 8],
               [14, 15],
               [8, 10],
               [8, 9],
               [8, 6],
               [10, 11],
               [10, 12],
               [10, 13],
               [7, 6],
               [14, 16],
               [18, 19],
               [18, 20],
               [18, 21],
               [18, 16],
               [17, 16]],
     'chains': [{'index': 0,
                 'residues': [{'atoms': [{'element': 'H',
                                          'index': 0,
                                          'name': 'H1'},
                                         {'element': 'C',
                                          'index': 1,
                                          'name': 'CH3'},
                                         {'element': 'H',
                                          'index': 2,
                                          'name': 'H2'},
                                         {'element': 'H',
                                          'index': 3,
                                          'name': 'H3'},
                                         {'element': 'C',
                                          'index': 4,
                                          'name': 'C'},
                                         {'element': 'O',
                                          'index': 5,
                                          'name': 'O'}],
                               'index': 0,
                               'resSeq': 1,
                               'name': 'ACE'},
                              {'atoms': [{'element': 'N',
                                          'index': 6,
                                          'name': 'N'},
                                         {'element': 'H',
                                          'index': 7,
                                          'name': 'H'},
                                         {'element': 'C',
                                          'index': 8,
                                          'name': 'CA'},
                                         {'element': 'H',
                                          'index': 9,
                                          'name': 'HA'},
                                         {'element': 'C',
                                          'index': 10,
                                          'name': 'CB'},
                                         {'element': 'H',
                                          'index': 11,
                                          'name': 'HB1'},
                                         {'element': 'H',
                                          'index': 12,
                                          'name': 'HB2'},
                                         {'element': 'H',
                                          'index': 13,
                                          'name': 'HB3'},
                                         {'element': 'C',
                                          'index': 14,
                                          'name': 'C'},
                                         {'element': 'O',
                                          'index': 15,
                                          'name': 'O'}],
                               'index': 1,
                               'resSeq': 2,
                               'name': 'ALA'},
                              {'atoms': [{'element': 'N',
                                          'index': 16,
                                          'name': 'N'},
                                         {'element': 'H',
                                          'index': 17,
                                          'name': 'H'},
                                         {'element': 'C',
                                          'index': 18,
                                          'name': 'C'},
                                         {'element': 'H',
                                          'index': 19,
                                          'name': 'H1'},
                                         {'element': 'H',
                                          'index': 20,
                                          'name': 'H2'},
                                         {'element': 'H',
                                          'index': 21,
                                          'name': 'H3'}],
                               'index': 2,
                               'resSeq': 3,
                               'name': 'NME'}]}]}

