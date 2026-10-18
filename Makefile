MODULE_TOPDIR = ../..

PGM = i.hyper.cloudremoval

ETCFILES = cloudlib solar_table

ifeq ($(DEBIAN_BUILD),1)
PREFIX   ?= /usr/local
BINDIR    = $(DESTDIR)$(PREFIX)/bin
LIBDIR    = $(DESTDIR)$(PREFIX)/lib/$(PGM)
PGM_PY   = $(PGM).py

$(PGM): $(PGM_PY)
	cp $< $@ && chmod 755 $@

install: $(PGM)
	install -d $(BINDIR) $(LIBDIR)
	install -m 755 $(PGM) $(BINDIR)/$(PGM)
	install -m 644 $(addsuffix .py,$(ETCFILES)) $(LIBDIR)/

clean:
	rm -f $(PGM)

.PHONY: install clean

else
include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Python.make

default: script
endif
