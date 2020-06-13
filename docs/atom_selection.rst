Atom Selection Reference
========================

Introduction
------------

MDTaj's :ref:`trajectory analysis <analysis>` functions use 0-based arrays
of "atom indices" to refer to subsets or groups of atoms in trajectories.
To generate these index arrays, MDTraj includes a powerful text-based atom
selection domain-specific language operating on the ``Topology``. The
following are all valid selection queries::

    top.select("water")
    top.select("resSeq 35")
    top.select("water and name O")
    top.select("mass 5.5 to 20")
    top.select("resname =~ 'C.*'")
    top.select("protein and (backbone or resname ALA)")

These queries return a numpy array of integers containing the indices of
the matching atoms. Equivalent python code for every selection expression
can be generated using ``Topology.select_expression``. ::

    >>> top.select_expression("water and name O")
    "[atom.index for atom in topology.atoms if (atom.residue.is_water and (atom.name == 'O'))]"

Reference
---------

MDTraj recognizes the following keywords. Each keyword maps directly to a
property on the MDTraj topology object's ``Atom``/``Residue``/``Chain`` tree.

=============    ========================   =========      ================================================================
Keyword          Synonyms                   Type           Description
-------------    ------------------------   ---------      ----------------------------------------------------------------
``all``          ``everything``             ``bool``       Matches everything
``none``         ``nothing``                ``bool``       Matches nothing
``backbone``     ``is_backbone``            ``bool``       Whether atom is in the backbone of a protein residue
``sidechain``    ``is_sidechain``           ``bool``       Whether atom is in the sidechain of a protein residue
``protein``      ``is_protein``             ``bool``       Whether atom is part of a protein residue
``water``        ``is_water``, ``waters``   ``bool``       Whether atom is part of a water residue
``name``                                    ``str``        Atom name
``index``                                   ``int``        Atom index (0-based)
``n_bonds``                                 ``int``        Number of bonds this atom participates in
``type``         ``element``, ``symbol``    ``str``        1 or 2-letter chemical symbols from the periodic table
``mass``                                    ``float``      Element atomic mass (daltons)
``residue``      ``resSeq``                 ``int``        Residue Sequence record (generally 1-based, but depends on topology)
``resid``        ``resi``                   ``int``        Residue index (0-based)
``resname``      ``resn``                   ``str``        Residue name
``rescode``      ``code``, ``resc```        ``str``        1-letter residue code
``chainid``                                 ``int``        Chain index (0-based)
=============    ========================   =========      ================================================================

Literals
~~~~~~~~

Integer, floating point, and string literals are also parsed. Both single-quoted,
strings, double-quoted strings, and bare words are also parsed as string
literals. ::

    # The following queries are equivalent
    top.select("symbol == O")
    top.select("symbol == 'O'")
    top.select('symbol == "O"')

Operators
~~~~~~~~~

Standard boolean operations (``and``, ``or``, and ``not``) as well as their
C-style aliases (``&&``, ``||``, ``!``) are supported. The expected logical
operators (``<``, ``<=``, ``==``, ``!=``, ``>=``, ``>``) are also available, as
along with their FORTRAN-style synonyms (``lt``, ``le``, ``eq``, ``ne``,
``ge``, ``gt``).

A regular-expression matching operator, ``=~``, is available. For example, to
match any of the names ``'C1'``, ``'C2'``, ``'C3'``, ``'C4'``, you can use the
following query. The regular expression syntax is just the `native python Regex
syntax <https://docs.python.org/3/library/re.html#regular-expression-syntax>`_ ::

    top.select("name =~ 'C[1-4]'")

An implicit equality relation is implied between adjacent expressions ::

    # The following queries are equivalent
    top.select("resid 35")
    top.select("resid == 35")

Range queries
~~~~~~~~~~~~~

Range queries are also supported. The range condition is an expression of
the form ``<expression> <low> to <high>``, which resolves to ``<low> <=
<expression> <= <high>``.  For example ::

    # The following queries are equivalent
    top.select("resid 10 to 30")
    top.select("(10 <= resid) and (resid <= 30)")


Implementation
--------------

MDTraj atom selection DSL lets users specify an expression which operates
on a single ``Atom`` and return a ``bool``, which is used subsequently as
the predicate for a ``filter`` expression. The expressions compile to
Python bytecode, and are then executed directly against the topology object
in the python VM.

This is done in two steps: first, query strings are parsed according to a
grammar defined using `PyParsing <http://pyparsing.wikispaces.com/>`_. The
parse tree is traversed, and used to construct an `abstract syntax tree
<https://docs.python.org/3/library/ast.html>`_ corresponding to the
equivalent Python atom selection expression (e.g.
``Topology.select_expression``).

.. vim: tw=75
