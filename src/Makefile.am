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
	svdcmp.h \
	update.cc \
	phase_table.cc \
	phase_enum.hh \
	model.hh \
	element.hh \
	constants.hh \
	model.cc \
	melt.cc \
	algorithm.hh \
	solid.cc \
	vapor.cc \
	pfitzer_sterner.cc \
	melt_model.cc \
	magma.cc \
	solid.hh \
	vapor.hh \
	pfitzer_sterner.hh \
	bm_melt.cc \
	melt_model.hh \
	bm_melt.hh \
	melt.hh \
	landau_solid.hh \
	landau_solid.cc \
	bragg_williams_solid.cc \
	bragg_williams_solid.hh \
	bm_solid.cc \
	bm_solid.hh

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

