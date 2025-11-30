## Copyright 2015-2016 Carnë Draug
## Copyright 2015-2016 Oliver Heimlich
## Copyright 2015-2019 Mike Miller
## Copyright 2024 John Donoghue
## Copyright 2025 Torsten Lilge
##
## Copying and distribution of this file, with or without modification,
## are permitted in any medium without royalty provided the copyright
## notice and this notice are preserved.  This file is offered as-is,
## without any warranty.

OCTAVE    ?= octave --silent
SED       := sed
MV        ?= mv -f
MAKEINFO  ?= makeinfo
# use --force in makeinfo since there are @refs to octave core functions
MAKEINFO_OPTIONS := --no-headers --no-split --no-validate \
												--set-customization-variable 'COPIABLE_LINKS 0'
# work out a possible help generator
ifeq ($(strip $(QHELPGENERATOR)),)
  ifneq ($(shell qhelpgenerator -v 2>/dev/null),)
    QHELPGENERATOR = qhelpgenerator
  else ifneq ($(shell qcollectiongenerator -v 2>/dev/null),)
    QHELPGENERATOR = qcollectiongenerator-qt5
  else ifneq ($(shell qcollectiongenerator -qt5 -v 2>/dev/null),)
    QHELPGENERATOR = qcollectiongenerator -qt5
  else
    QHELPGENERATOR = true
  endif
endif

DESCRIPTION := DESCRIPTION
PACKAGE := $(shell $(SED) -n -e 's/^[Nn]ame: *\(\w\+\)/\1/p' $(DESCRIPTION))
VERSION := $(shell $(SED) -n -e 's/^[Vv]ersion: *\(\w\+\)/\1/p' $(DESCRIPTION))
DATE    := $(shell $(SED) -n -e 's/^[Dd]ate: *\(\w\+\)/\1/p' $(DESCRIPTION))
DEPENDS := $(shell $(SED) -n -e 's/^[Dd]epends[^,]*, *\(.*\)/\1/p' $(DESCRIPTION) | $(SED) 's/ *([^()]*)//g; s/ *, */ /g')

TARGET_DIR      := target
RELEASE_DIR     := $(TARGET_DIR)/$(PACKAGE)-$(VERSION)
RELEASE_TARBALL := $(TARGET_DIR)/$(PACKAGE)-$(VERSION).tar.gz
DOCS_HTML_DIR   := docs
DOCS_DIR        := doc
DOCS_DEV_DIR    := devel/doc

SC_SUBMOD       := slicot-reference
SRC             := src
SC              := $(SRC)/slicot
SC_SRC          := src
SC_LAPACK       := $(SC_SRC)/lapack_aux
SC_DOC          := doc/SLICOT

M_SOURCES       := $(wildcard inst/*.m)
CC_SOURCES      := $(wildcard src/*.cc)
TEXI_TMP        := $(DOCS_DEV_DIR)/functions.texi $(DOCS_DEV_DIR)/version.texi
LOGO_DEV        := $(DOCS_DEV_DIR)/$(PACKAGE)-logo.svg
DOCS_SOURCES    := $(DOCS_DEV_DIR)/$(PACKAGE).texi $(DOCS_DEV_DIR)/copying.texi \
                   $(LOGO_DEV) $(TEXI_TMP)
PKG_ADD         := $(shell grep -sPho '(?<=(//|\#\#) PKG_ADD: ).*' $(CC_SOURCES) $(M_SOURCES))

DOCS_PDF        := $(DOCS_DIR)/$(PACKAGE).pdf
DOCS_QCH        := $(DOCS_DIR)/$(PACKAGE).qch
DOCS_LOGO       := $(DOCS_DIR)/$(PACKAGE).svg

.PHONY: help dist docs-html docs release install all check run clean

help:
	@echo " "
	@echo "Targets:"
	@echo " "
	@echo "   dist      - Create $(RELEASE_TARBALL) for release"
	@echo "   docs-html - Create html documentation in $(DOCS_DIR), assumes"
	@echo "               current version of $(PACKAGE) is installed;"
	@echo "   docs        Create the pdf manual and the Qt help file"
	@echo "               in $(DOC_DIR)"
	@echo "   qch       - Build Qt help file"
	@echo "   pdf       - Build pdf manual"
	@echo "   release   - Create dist, install and create docs,"
	@echo "               shows sha256sum of tarball"
	@echo
	@echo "   install   - Install the package in GNU Octave"
	@echo "   all       - Build all oct files"
	@echo "   check     - Execute package tests (with install)"
	@echo
	@echo "   clean     - Remove releases, doc and oct files"
	@echo "   distclean - Remove releases, oct files and compiled libraries"
	@echo " "

%.tar.gz: %
	@echo "Create $@ ..."
	@tar -c -f - --posix -C "$(TARGET_DIR)/" "$(notdir $<)" | gzip -9n > "$@"

# Create the directory structure of the destributed archive file.
# For this, archive the current git repo with some exceptions given
# in .gitattributes with export-ignore and untar it. Then, copy the
# slicot routines together with README and LICENSE file into a
# subdirectory. Finally move  TG04BX.f, a version of TB04BX.f extended
# for descriptor systems into that directory.
$(RELEASE_DIR): .git/index
	@echo "Creating package dist directory $@ ..."
	@-$(RM) -r $@
	@mkdir -p $@
	@echo "  git archive ..."
	@git archive -o $@/tmp.tar HEAD
	@cd $@ && tar -xf tmp.tar && $(RM) tmp.tar
	@echo "  copy doc files ..."
	@cp $(DOCS_PDF) $@/$(DOCS_DIR)/
	@cp $(DOCS_QCH) $@/$(DOCS_DIR)/
	@cp $(DOCS_LOGO) $@/$(DOCS_DIR)/
	@echo "  copy slicot files ..."
	@mkdir -p $@/$(SC)/$(SC_LAPACK)
	@cp -t $@/$(SC)/$(SC_SRC) $(SC_SUBMOD)/$(SC_SRC)/*.f
	@cp -t $@/$(SC)/$(SC_LAPACK) $(SC_SUBMOD)/$(SC_LAPACK)/*.f
	@cp $(SC_SUBMOD)/LICENSE   $@/$(SC_SRC)/../
	@cp $(SC_SUBMOD)/README.md $@/$(SC_SRC)/../README-SLICOT.md
	@cp $(SC_SUBMOD)/LICENSE   $@/$(SC_DOC)/
	@cp $(SC_SUBMOD)/README.md $@/$(SC_DOC)/README-SLICOT.md
	@echo "  bootstrap ..."
	@cd $@/$(SRC) && ./bootstrap && $(RM) -r "autom4te.cache"
	@chmod -R a+rX,u+w,go-w "$@"

docs-html:
	@echo "Updating HTML documentation ... "
	@cd $(DOCS_HTML_DIR) && $(OCTAVE) \
													--eval "pkg load pkg-octave-doc; " \
													--eval "pkg load $(PACKAGE);" \
													--eval 'package_texi2html ("${PACKAGE}");'

$(DOCS_DEV_DIR)/version.texi: $(DESCRIPTION)
	@echo Generating $@ ...
	@echo "@c autogenerated from Makefile" > $@
	@echo "@set VERSION $(VERSION)" >> $@
	@echo "@set PACKAGE $(PACKAGE)" >> $@
	@echo "@set DATE $(DATE)" >> $@

# Collect all texinfo tests in one file. For being able to use
# makeinfo --pdf, --html and pdf-octave-doc and always have nice
# refs and links from @seealso commands in all output formats.
$(DOCS_DEV_DIR)/functions.texi: .git/index
	@echo Generating $@ \(collect all texinfo texts\) ...
	@$(DOCS_DEV_DIR)/mkfuncdocs.py --src-dir=inst/ --src-dir=src/ INDEX > $@
	@$(DOCS_DEV_DIR)/fix_seealso.pl $@

$(DOCS_LOGO): $(LOGO_DEV)
	@echo Copying $@ ...
	cp $(LOGO_DEV) $@

$(DOCS_PDF): $(DOCS_SOURCES)
	@echo Generating $@ ...
	@cd $(DOCS_DEV_DIR) && $(MAKEINFO) --pdf $(MAKEINFO_OPTIONS) $(PACKAGE).texi > /dev/null 2>&1
	@cd $(DOCS_DEV_DIR) && $(RM) *.out *.log *.idx *.ilg *.ind *.toc *.cp *.cps *.aux *.fn *.fns *.out
	@$(MV) $(DOCS_DEV_DIR)/$(PACKAGE).pdf $(DOCS_PDF)

$(DOCS_QCH): $(DOCS_SOURCES)
	@echo Generating $@ ...
	@cd $(DOCS_DEV_DIR) && $(MAKEINFO) --html --css-ref=$(PACKAGE).css $(MAKEINFO_OPTIONS) $(PACKAGE).texi
ifeq ($(QHELPGENERATOR),true)
	$(warning No QHELPGENERATOR... skipping QT doc build)
else
	@cd $(DOCS_DEV_DIR) && ./mkqhcp.py $(PACKAGE)\
									&& $(QHELPGENERATOR) -s $(PACKAGE).qhcp -o $(PACKAGE).qhc > /dev/null 2>&1
	@cd $(DOCS_DEV_DIR) && $(RM) $(PACKAGE).qhcp $(PACKAGE).qhp $(PACKAGE).qhc $(PACKAGE).html
endif
	@$(MV) $(DOCS_DEV_DIR)/$(PACKAGE).qch $(DOCS_QCH)

qch: $(DOCS_QCH)

pdf: $(DOCS_PDF)

docs: qch pdf $(DOCS_LOGO)

dist: docs $(RELEASE_TARBALL)

release: dist install
	sha256sum $(RELEASE_TARBALL)
	@echo ''
	@echo '* Execute: git tag "control-${VERSION}"'
	@echo "* Push release branch to https://github.com/gnu-octave/pkg-control"
	@echo "* Push new tag to https://github.com/gnu-octave/pkg-control"
	@echo "* Create a new release on https://github.com/gnu-octave/pkg-control"
	@echo "* Update control.yaml in https://github.com/gnu-octave/packages"
	@echo "  and create/merge a pull request for the changes in this file"

install: dist $(RELEASE_TARBALL)
	@echo 'Installing package "${RELEASE_TARBALL}" locally ...'
	@$(OCTAVE) --eval 'pkg ("install", "${RELEASE_TARBALL}")'

all: $(CC_SOURCES)
	$(MAKE) -C $(SRC) all

check: install
	$(OCTAVE) --path "inst/" --path "src/" \
	  --eval 'pkg test control'

clean:
	$(RM) -r $(TARGET_DIR)
	$(MAKE) -C $(SRC) clean
	$(RM) $(DOCS_PDF) $(DOCS_QCH)
	$(RM) $(TEXI_TMP)

distclean: clean
	$(MAKE) -C $(SRC) distclean
