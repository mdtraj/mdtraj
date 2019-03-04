# Copied from the yt_project, commit e8fb57e
# yt/doc/extensions/notebook_sphinxext.py
#  https://bitbucket.org/yt_analysis/yt/src/e8fb57e66ca42e26052dadf054a5c782740abec9/doc/extensions/notebook_sphinxext.py?at=yt

# Almost completely re-written by Matthew Harrigan to use nbconvert v4

from __future__ import print_function

import os
import shutil

from docutils import nodes
from docutils.parsers.rst import directives, Directive
import nbformat
from nbconvert import HTMLExporter, PythonExporter


def _read(wd, name):
    with open("{}/{}.ipynb".format(wd, name)) as f:
        notebook = nbformat.read(f, as_version=4)
    return notebook


def export_html(wd, name):
    nb = _read(wd, name)

    config = {
        'Exporter': {'template_file': 'embed',
                     'template_path': ['./sphinxext/']},
        'ExecutePreprocessor': {'enabled': True},
        'ExtractOutputPreprocessor': {'enabled': True},
        'CSSHTMLHeaderPreprocessor': {'enabled': True}
    }

    exporter = HTMLExporter(config)

    try:
        body, resources = exporter.from_notebook_node(nb)

        for fn, data in resources['outputs'].items():
            with open("{}/{}".format(wd, fn), 'wb') as f:
                f.write(data)
        return body
    except Exception as e:
        return str(e)


def export_python(wd, name):
    nb = _read(wd, name)
    exporter = PythonExporter()
    body, resources = exporter.from_notebook_node(nb)
    with open("{}/{}.py".format(wd, name), 'w') as f:
        f.write(body)


class NotebookDirective(Directive):
    """Insert an evaluated notebook into a document
    """
    required_arguments = 1
    optional_arguments = 1
    option_spec = {'skip_exceptions': directives.flag}
    final_argument_whitespace = True

    def run(self):

        # check if raw html is supported
        if not self.state.document.settings.raw_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)

        # get path to notebook
        nb_rel_path = self.arguments[0]
        nb_abs_path = "{}/../{}".format(setup.confdir, nb_rel_path)
        nb_abs_path = os.path.abspath(nb_abs_path)
        nb_name = os.path.basename(nb_rel_path).split(".")[0]
        dest_dir = "{}/{}/{}".format(
            setup.app.builder.outdir,
            os.path.dirname(nb_rel_path),
            nb_name)
        fmt = {'wd': dest_dir, 'name': nb_name}

        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        shutil.copyfile(nb_abs_path, "{wd}/{name}.ipynb".format(**fmt))

        # TODO: Actually save evaluated notebook
        shutil.copyfile(nb_abs_path, "{wd}/{name}_eval.ipynb".format(**fmt))

        html = export_html(**fmt)
        export_python(**fmt)

        # Create link to notebook and script files
        link_rst = "({uneval}; {eval}; {py})".format(
            uneval=formatted_link("{wd}/{name}.ipynb".format(**fmt)),
            eval=formatted_link("{wd}/{name}_eval.ipynb".format(**fmt)),
            py=formatted_link("{wd}/{name}.py".format(**fmt)),
        )

        rst_file = self.state_machine.document.attributes['source']
        self.state_machine.insert_input([link_rst], rst_file)

        # create notebook node
        attributes = {'format': 'html', 'source': 'nb_path'}
        nb_node = notebook_node('', html, **attributes)
        nb_node.source, nb_node.line = self.state_machine \
            .get_source_and_line(self.lineno)

        # add dependency
        self.state.document.settings.record_dependencies.add(nb_abs_path)

        return [nb_node]


class notebook_node(nodes.raw):
    pass


def formatted_link(path):
    return "`%s <%s>`__" % (os.path.basename(path), path)


def visit_notebook_node(self, node):
    self.visit_raw(node)


def depart_notebook_node(self, node):
    self.depart_raw(node)


def setup(app):
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    app.add_node(notebook_node,
                 html=(visit_notebook_node, depart_notebook_node))

    app.add_directive('notebook', NotebookDirective)
