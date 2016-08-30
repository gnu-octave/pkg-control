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
HTML_DIR        := $(TARGET_DIR)/$(PACKAGE)-html
HTML_TARBALL    := $(TARGET_DIR)/$(PACKAGE)-html.tar.gz

M_SOURCES   := $(wildcard inst/*.m)
CC_SOURCES  := $(wildcard src/*.cc)
PKG_ADD     := $(shell grep -sPho '(?<=(//|\#\#) PKG_ADD: ).*' $(CC_SOURCES) $(M_SOURCES))

OCTAVE ?= octave --no-window-system --silent

.PHONY: help dist html release install all check run clean

help:
	@echo "Targets:"
	@echo "   dist    - Create $(RELEASE_TARBALL) for release"
	@echo "   html    - Create $(HTML_TARBALL) for release"
	@echo "   release - Create both of the above and show md5sums"
	@echo
	@echo "   install - Install the package in GNU Octave"
	@echo "   all     - Build all oct files"
	@echo "   check   - Execute package tests (w/o install)"
	@echo "   run     - Run Octave with development in PATH (no install)"
	@echo
	@echo "   clean   - Remove releases, html documentation, and oct files"

%.tar.gz: %
	tar -c -f - --posix -C "$(TARGET_DIR)/" "$(notdir $<)" | gzip -9n > "$@"

$(RELEASE_DIR): .hg/dirstate
	@echo "Creating package version $(VERSION) release ..."
	-$(RM) -r "$@"
	hg archive \
	  --exclude ".hg*" --exclude "Makefile" --exclude "devel" \
	  --type files "$@"
	chmod -R a+rX,u+w,go-w "$@"

$(HTML_DIR): install
	@echo "Generating HTML documentation. This may take a while ..."
	-$(RM) -r "$@"
	$(OCTAVE) --no-window-system --silent \
	  --eval "pkg load generate_html; " \
	  --eval "pkg load $(PACKAGE);" \
	  --eval 'generate_package_html ("${PACKAGE}", "$@", "octave-forge");'
	chmod -R a+rX,u+w,go-w $@

dist: $(RELEASE_TARBALL)
html: $(HTML_TARBALL)

release: dist html
	md5sum $(RELEASE_TARBALL) $(HTML_TARBALL)
	@echo "Upload @ https://sourceforge.net/p/octave/package-releases/new/"
	@echo 'Execute: hg tag "release-${VERSION}"'

## Note that in development versions this target may fail if we are dependent
## on unreleased versions.  This is by design, to force possible developers
## to set this up by hand (either using the "-nodeps" option" or changing the
## dependencies on DESCRIPTION.
install: $(RELEASE_TARBALL)
	@echo "Installing package locally ..."
	$(OCTAVE) --eval 'pkg ("install", "${RELEASE_TARBALL}")'

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

distclean: clean
	$(MAKE) -C src/ distclean
