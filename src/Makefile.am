## Process this file with automake to produce Makefile.in

## Created by Anjuta


uidir = $(pkgdatadir)/ui
ui_DATA = melt.ui


AM_CPPFLAGS = \
	-DPACKAGE_LOCALE_DIR=\""$(localedir)"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir). \
	$(ds__-0.1_CFLAGS). \
	$(gsl_CFLAGS). \
	$(ds__-0.1_CFLAGS). \
	$(ds__-0.1_CFLAGS). \
	$(libkgb-0.1_CFLAGS). \
	$(LIBKGB_CFLAGS). \
	$(libkgb-0.1_CFLAGS)"\" \
	-DPACKAGE_DATA_DIR=\""$(pkgdatadir)"\" \
	$(MELT_CFLAGS)

AM_CFLAGS =\
	 -Wall\
	 -g

bin_PROGRAMS = melt

melt_SOURCES = \
	main.cc \
	phase.hh \
	Gf.cc \
	phase.cc \
	svdcmp.h

melt_LDFLAGS = 

melt_LDADD = $(MELT_LIBS) \
	$(libkgb-0.1_LIBS) \
	$(LIBKGB_LIBS) \
	$(ds__-0.1_LIBS) \
	$(gsl_LIBS)




EXTRA_DIST = $(ui_DATA)

# Remove ui directory on uninstall
uninstall-local:
	-rm -r $(uidir)
	-rm -r $(pkgdatadir)

