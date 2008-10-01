sinclude ../../Makeconf

PKG_FILES = COPYING DESCRIPTION INDEX $(wildcard inst/*) doc/control.pdf doc/control.texi doc/control.txi
SUBDIRS = doc/

.PHONY: $(SUBDIRS)

pre-pkg::
	@for _dir in $(SUBDIRS); do \
	  $(MAKE) -C $$_dir all; \
	done

clean:
	@for _dir in $(SUBDIRS); do \
	  $(MAKE) -C $$_dir $(MAKECMDGOALS); \
	done
