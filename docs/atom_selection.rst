Atom Selection Language
=======================

Introduction
------------

MDTraj includes a powerful text-based atom selection / query language. The
following are all valid selection queries::

    traj.select("water and name O")
    traj.select("mass 5.5 to 20")
    traj.select("protein and (backbone or resname ALA)")

Valid boolean operations are ``and``, ``or``, and ``not`` as well as their
C-style aliases (``&&``, ``||``, ``!``). Comparison operators are also
allowed (defaults to ``==`` if not specified). Consider the following
equivalencies::

    # The following two queries are equivalent
    traj.select("resid > 200")
    traj.select("resid gt 200")

    # These are also equivalent
    traj.select("resname ALA")
    traj.select("resname == 'ALA'")


Reference
---------
MDTraj recognizes the following keywords. "Property Name" is the internal
name used by MDTraj and its API.

Per-Atom Fields
~~~~~~~~~~~~~~~

.. atom_select_table:: AtomUnaryOperand

.. atom_select_table:: AtomBinaryOperand

.. atom_select_table:: ElementBinaryOperand


Per-Residue Fields
~~~~~~~~~~~~~~~~~~

.. atom_select_table:: ResidueUnaryOperand

.. atom_select_table:: ResidueBinaryOperand


Advanced Usage
--------------

The atom selection language is very powerful and should be sufficient
for most users. For users who would benefit from a more direct interaction
with the MDTraj API, the atom selection language can be used to generate
syntactically valid python expressions that can be embedded in scripts or
applications built on MDTraj::

    traj.select_expression("water and (not name H)", top_name="traj.topology")
    "[a.index for a in traj.topology.atoms if (a.residue.is_water and (not a.name == 'H'))]"
