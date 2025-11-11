## Copyright 2015-2016 Carnë Draug
## Copyright 2015-2016 Oliver Heimlich
##
## Copying and distribution of this file, with or without modification,
## are permitted in any medium without royalty provided the copyright
## notice and this notice are preserved.  This file is offered as-is,
## without any warranty.

PACKAGE := $(shell grep "^Name: " DESCRIPTION | cut -f2 -d" ")
VERSION := $(shell grep "^Version: " DESCRIPTION | cut -f2 -d" ")

TARGET_DIR      := target
RELEASE_DIR     := $(TARGET_DIR)/$(PACKAGE)-$(VERSION)
RELEASE_TARBALL := $(TARGET_DIR)/$(PACKAGE)-$(VERSION).tar.gz
DOCS_DIR        := docs
DOC_DIR         := doc

SC_SMOD := slicot-reference
SC_SRC  := src/slicot/src
SC_DOC  := doc/SLICOT

M_SOURCES   := $(wildcard inst/*.m)
CC_SOURCES  := $(wildcard src/*.cc)

PKG_ADD     := $(shell grep -sPho '(?<=(//|\#\#) PKG_ADD: ).*' $(CC_SOURCES) $(M_SOURCES))

OCTAVE ?= octave --silent

.PHONY: help dist docs-html doc release install all check run clean

test:
	@echo $(SLICOT_SOURCES)

help:
	@echo " "
	@echo "Targets:"
	@echo " "
	@echo "   dist      - Create $(RELEASE_TARBALL) for release"
	@echo "   docs-html - Create html documentation in $(DOCS_DIR), assumes"
	@echo "               current version of $(PACKAGE) is installed;"
	@echo "   doc         Create the pdf manual and the Qt help file"
	@echo "               in $(DOC_DIR)"
	@echo "   release   - Create dist, install and create docs,"
	@echo "               shows sha256sum of tarball"
	@echo
	@echo "   install   - Install the package in GNU Octave"
	@echo "   all       - Build all oct files"
	@echo "   check     - Execute package tests (w/o install)"
	@echo "   run       - Run Octave with development in PATH (no install)"
	@echo
	@echo "   clean     - Remove releases and oct files"
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
	@echo "  copy packinfo ..."
	@echo "  copy pdf and qch file ..."
	@cp $(DOC_DIR)/$(PACKAGE).pdf $@/$(DOC_DIR)/
	@cp $(DOC_DIR)/$(PACKAGE).qch $@/$(DOC_DIR)/
	@echo "  copy files from slicot ..."
	@mkdir -p $@/$(SC_SRC)
	@cp -t $@/$(SC_SRC) $(SC_SMOD)/src/*.f
	@cp -t $@/$(SC_SRC) $(SC_SMOD)/src_aux/*.f
	@cp $(SC_SMOD)/LICENSE   $@/$(SC_SRC)/../
	@cp $(SC_SMOD)/README.md $@/$(SC_SRC)/../README-SLICOT.md
	@cp $(SC_SMOD)/LICENSE   $@/$(SC_DOC)/
	@cp $(SC_SMOD)/README.md $@/$(SC_DOC)/README-SLICOT.md
	@echo "  bootstrap ..."
	@cd $@/src && ./bootstrap && $(RM) -r "autom4te.cache"
	@chmod -R a+rX,u+w,go-w "$@"

docs-html:
	@echo "Updating HTML documentation. This may take a while ..."
	@cd "$(DOCS_DIR)" && $(OCTAVE) \
	                      --eval "pkg load pkg-octave-doc; " \
	                      --eval "pkg load $(PACKAGE);" \
	                      --eval 'package_texi2html ("${PACKAGE}");'

doc:
	@echo "Creating pdf and qch file in $(DOC_DIR) ..."
	@cd "$(DOC_DIR)" && make

dist: doc $(RELEASE_TARBALL)

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
	$(MAKE) -C src/ all

check: all
	$(OCTAVE) --path "inst/" --path "src/" \
	  --eval '${PKG_ADD}' \
	  --eval 'runtests ("inst"); runtests ("src");'

run: all
	$(OCTAVE) --persist --path "inst/" --path "src/" \
	  --eval '${PKG_ADD}'

clean:
	$(RM) -r $(TARGET_DIR)
	$(MAKE) -C src/ clean
	@cd "$(DOC_DIR)" && make clean


distclean: clean
	$(MAKE) -C src/ distclean
