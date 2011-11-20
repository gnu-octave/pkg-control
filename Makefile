sinclude ../../Makeconf

PKG_FILES = COPYING DESCRIPTION INDEX NEWS $(wildcard inst/*)
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
