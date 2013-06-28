# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = ../../mdtraj-docs

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

#.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext

all: html

view:
	cd .. && python setup.py install
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)"

github: html
	cd $(BUILDDIR); \
	git commit -am "new docs built"; \
	git push origin gh-pages


clean:
	mv $(BUILDDIR)/.git tmp.git
	mv $(BUILDDIR)/.nojekyll tmp.nojekyll
	rm -rf $(BUILDDIR)/*
	mv tmp.git $(BUILDDIR)/.git
	mv tmp.nojekyll $(BUILDDIR)/.nojekyll 
	
# dirhtml:
# 	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
# 	@echo
# 	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

# singlehtml:
# 	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
# 	@echo
# 	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

# pickle:
# 	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
# 	@echo
# 	@echo "Build finished; now you can process the pickle files."

# json:
# 	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
# 	@echo
# 	@echo "Build finished; now you can process the JSON files."

# htmlhelp:
# 	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
# 	@echo
# 	@echo "Build finished; now you can run HTML Help Workshop with the" \
# 	      ".hhp project file in $(BUILDDIR)/htmlhelp."

# qthelp:
# 	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(BUILDDIR)/qthelp
# 	@echo
# 	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
# 	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
# 	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/MDTraj.qhcp"
# 	@echo "To view the help file:"
# 	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/MDTraj.qhc"

# devhelp:
# 	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(BUILDDIR)/devhelp
# 	@echo
# 	@echo "Build finished."
# 	@echo "To view the help file:"
# 	@echo "# mkdir -p $$HOME/.local/share/devhelp/MDTraj"
# 	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/MDTraj"
# 	@echo "# devhelp"

# epub:
# 	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
# 	@echo
# 	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

# latex:
# 	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
# 	@echo
# 	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
# 	@echo "Run \`make' in that directory to run these through (pdf)latex" \
# 	      "(use \`make latexpdf' here to do that automatically)."

# latexpdf:
# 	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
# 	@echo "Running LaTeX files through pdflatex..."
# 	$(MAKE) -C $(BUILDDIR)/latex all-pdf
# 	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

# text:
# 	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(BUILDDIR)/text
# 	@echo
# 	@echo "Build finished. The text files are in $(BUILDDIR)/text."

# man:
# 	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(BUILDDIR)/man
# 	@echo
# 	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

# texinfo:
# 	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
# 	@echo
# 	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
# 	@echo "Run \`make' in that directory to run these through makeinfo" \
# 	      "(use \`make info' here to do that automatically)."

# info:
# 	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
# 	@echo "Running Texinfo files through makeinfo..."
# 	make -C $(BUILDDIR)/texinfo info
# 	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

# gettext:
# 	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
# 	@echo
# 	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

# changes:
# 	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
# 	@echo
# 	@echo "The overview file is in $(BUILDDIR)/changes."

# linkcheck:
# 	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
# 	@echo
# 	@echo "Link check complete; look for any errors in the above output " \
# 	      "or in $(BUILDDIR)/linkcheck/output.txt."

# doctest:
# 	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
# 	@echo "Testing of doctests in the sources finished, look at the " \
# 	      "results in $(BUILDDIR)/doctest/output.txt."
