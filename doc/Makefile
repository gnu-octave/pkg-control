sinclude ../../../Makeconf

# Fill in the variables as it makes testing the package manager easier
ifeq ($(MKDOC),)
MKDOC = ../../../admin/mkdoc
MKTEXI = ../../../admin/mktexi
MAKEINFO = makeinfo --no-split
TEXI2PDF = texi2ppf --clean
DVIPS = dvips
LN_S = ln -s
endif

INFODOC = control.info
PSDOC = $(patsubst %.info,%.ps,$(INFODOC))
PDFDOC = $(patsubst %.info,%.pdf,$(INFODOC))
HTMLDOC = $(patsubst %.info,%.html,$(INFODOC))
TEXIDOC = $(patsubst %.info,%.texi,$(INFODOC))
DOCS = $(INFODOC) $(PDFDOC)
DOCSTRINGS = DOCSTRINGS
INDEX = ../INDEX
TMPDELETES = *.log *.dvi $(DOCSTRINGS) $(TEXIDOC) *~
DELETES = $(TMPDELETES) *.ps *.texi *.info $(DOCS) *.html control/ html/

all : $(PDFDOC) $(HTMLDOC) ../inst/doc.info

../inst/doc.info : $(INFODOC)
	cp -f $(INFODOC) ../inst/doc.info

%.dvi : %.texi
	@if test "x$(TEXI2DVI)" != "x"; then \
	  echo "Making dvi $@"; \
	  TEXINPUTS="./:../../..:$(TEXINPUTS):"; \
          export TEXINPUTS; \
	  $(TEXI2DVI) $< ; \
	fi

%.ps : %.dvi 
	@if test "x$(TEXI2DVI)" != "x" && test "x$(DVIPS)" != "x"; then \
	  echo "Making postscript $@"; \
	  $(DVIPS) -o $@ $< ; \
	fi

ifeq (,$(TEXI2PDF))
%.pdf : %.dvi
	@if test "x$(TEXI2DVI)" != "x" && test "x$(DVIPDF)" != "x"; then \
	  echo "Making pdf $@"; \
	  $(DVIPDF) $< ; \
	fi
else
%.pdf : %.texi
	@if test "x$(TEXI2PDF)" != "x"; then \
	  echo "Making pdf $@"; \
	  TEXINPUTS="./:../../..:$(TEXINPUTS):"; \
          export TEXINPUTS; \
	  $(TEXI2PDF) $< ; \
	fi
endif

%.info : %.texi
	@if test "x$(MAKEINFO)" != "x"; then \
	  echo "Making info $@"; \
	  $(MAKEINFO) -I./ -I../../../ $< ; \
	fi

# This build target is broken as the refcard card is written with an old
# version of lout
csrefcard.ps : csrefcard.lt OCSTreport
	lout csrefcard.lt > csrefcard.ps

# Need a stupid copy of the TOC for older texi2html versions
# Newer texi2html place documentation in a sub-directory
%.html : %.texi
	@if test "x$(TEXI2HTML)" != "x"; then \
	  echo "Making html $@"; \
	  $(TEXI2HTML) -I . -I ../../.. -iftex --subdir=./ -expandinfo $< ; \
	  if test ! -e "$(@:.html=_toc.html)"; then \
	    if test ! -e "control/$@"; then \
	      $(INSTALL_DATA) control/$(@:.html=_toc.html) control/$@ ; \
	    fi; \
	    $(LN_S) $@ control/index.html; \
	    mv control html; \
	    $(INSTALL_DATA) html; \
	  else \
	    if test ! -e "$@"; then \
	      $(INSTALL_DATA) $(@:.html=_toc.html) $@ ; \
	    fi; \
	    if [ ! -e "html/" ]; then \
	      mkdir html; \
	    fi; \
	    $(INSTALL_DATA) *.html html; \
	    $(LN_S) $@ html/index.html; \
	  fi \
	fi

.PRECIOUS: %.texi
%.texi : %.txi
	@echo "Making texinfo $@"; \
	$(RM) -f $(DOCSTRINGS); \
	$(MKDOC) ../ > $(DOCSTRINGS); \
	$(MKTEXI) $< $(DOCSTRINGS) $(INDEX) > $@ ; \
	$(RM) -f $(DOCSTRINGS);

clean:
	@echo "Cleaning..."; \
	$(RM) -fr $(DELETES)

dist: all

count:
	wc *.txi

