__author__ = 'harrigan'

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives.tables import Table
from mdtraj.core import selection, topology

from docutils.nodes import tbody, entry, paragraph, row, table, tgroup, colspec, \
    thead


def rev_dic(dic, value):
    """Yield keys corresponding to a given value."""
    for k in dic:
        if dic[k] == value:
            yield k


def fields_str(aliases, field_name):
    """Make a list of aliases for a given field."""
    return ", ".join(rev_dic(aliases, field_name))


def description(top_class, field_name):
    """Get the first line of the docstring for each property."""
    try:
        desc = getattr(top_class, field_name).__doc__.splitlines()[0]
    except Exception as e:
        desc = "Couldn't get description: {}".format(e)
    return desc


def make_table_header():
    """Make a header row."""
    return thead(
        "", row(
            "",
            entry("", paragraph(text="Keywords")),
            entry("", paragraph(text="Property Name")),
            entry("", paragraph(text="Description"))
        )
    )


def make_table_body(aliases, top_class):
    """Make one entry for each field.

    Parameters
    ----------
    aliases : dict
        Dictionary of alias_keyword -> main_keyword
    """

    # Unique values define the actual different properties
    field_names = sorted(set(aliases.values()))

    my_tbody = tbody()

    # Add header

    # Add rows
    my_tbody += [
        row("",
            entry("", paragraph(text=fields_str(aliases, field_name))),
            entry("", paragraph(text=field_name)),
            entry("", paragraph(text=description(top_class, field_name)))
        )
        for field_name in field_names
    ]
    return my_tbody


def make_table(aliases, top_class):
    """Make a table node."""
    my_table = table(
        "",
        tgroup(
            "",
            colspec(colwidth=1),
            colspec(colwidth=1),
            colspec(colwidth=1),
            make_table_header(),
            make_table_body(aliases, top_class),
            cols=3
        )
    )
    return my_table


class AtomSelectDirective(Table):
    has_content = True
    final_argument_whitespace = True
    required_arguments = 1
    optional_arguments = 1

    def run(self):
        class_name = self.arguments[0]
        if len(self.arguments) == 2:
            table_desc = self.arguments[1]

        # Get class from classname
        cls = getattr(selection, class_name)

        # Get keywords from class (specified as input argument)
        aliases = cls.keyword_aliases

        # Get topology class (Residue, Atom, ...)
        top_name = cls.get_top_name()
        top_class = getattr(topology, top_name)

        return [make_table(aliases, top_class)]


def setup(app):
    app.add_directive('atom_select_table', AtomSelectDirective)