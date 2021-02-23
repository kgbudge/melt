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
	svdcmp.h \
	phase_table.cc \
	phase_enum.hh \
	element.hh \
	constants.hh \
	algorithm.hh \
	solid.cc \
	vapor.cc \
	pfitzer_sterner.cc \
	magma.cc \
	solid.hh \
	vapor.hh \
	pfitzer_sterner.hh \
	bm_melt.cc \
	bm_melt.hh \
	landau_solid.hh \
	landau_solid.cc \
	bragg_williams_solid.cc \
	bragg_williams_solid.hh \
	bm_solid.cc \
	bm_solid.hh \
	update.hh \
	element.cc \
	State.hh \
	D2.hh \
	D2.cc \
	on_new.cc \
	on_open.cc \
	on_quit.cc \
	on_save.cc \
	on_save_as.cc \
	on_activate.cc \
	on_activate_P.cc \
	on_activate_depth.cc \
	on_normalize.cc \
	update.cc \
	Model.cc \
	Model.hh \
	gui.hh \
	do_save.cc \
	on_by_toggle.cc \
	get_molar_composition.cc \
	get_composition.cc \
	iugs_classify.cc \
	classify.hh \
	classify_components.cc \
	State.cc \
	D1.cc \
	D1.hh \
	Phase.hh \
	Melt_Model.cc \
	Melt_Model.hh \
	Melt.hh \
	Melt.cc \
	Melt_Model__dGf.cc \
	State__melt.cc \
	Melt_Model__minimize_Gf.cc \
	State__update.cc \
	Melt_Model__minimize_trial_set_.cc \
	Melt_Model__compute_current_melt_composition_.cc \
	Melt_Model__Gfmelt.cc \
	Melt_Model__compute_melt_endmember_composition_.cc \
	Melt_Model__Gfm.cc \
	Melt_Model__compute_current_solid_composition_.cc \
	Melt_Model__Gf.cc \
	Melt_Model__update_current_state.cc \
	Melt_Model__update_solid_state_.cc \
	Melt_Model__construct_reaction_.cc \
	Melt_Model__dGfm.cc \
	Melt_Model__calculate_extent_.cc \
	Melt_Model__print.cc

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

