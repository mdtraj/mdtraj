__author__ = 'harrigan'

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives.tables import Table


class AtomSelectDirective(Table):
    has_content = True

    def run(self):
        table = nodes.table()
        tgroup = nodes.tgroup("", nodes.colspec(colwidth=1),
                              nodes.colspec(colwidth=1), cols=2)

        tbody = nodes.tbody()
        entry1 = nodes.entry("", nodes.paragraph(text="hi 1"))
        entry2 = nodes.entry("", nodes.paragraph(text="hi 2"))

        row1 = nodes.row("", entry1, entry2)

        tbody += row1
        tgroup += tbody
        table += tgroup
        return [table]


def setup(app):
    app.add_directive('atom_select_table', AtomSelectDirective)