/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * on_save_as.cc
 * Copyright (C) 2015 Kent G. Budge <kgb@kgbudge.com>
 * 
 * melt is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * melt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#if 0
#include <gtkmm.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "kgb/tostring.h"
#include "kgb/Assert.h"

#include "config.h"


#ifdef ENABLE_NLS
#  include <libintl.h> 
#endif

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"


/* For testing propose use the local (not installed) ui file */
#ifndef NDEBUG
#define UI_FILE PACKAGE_DATA_DIR"/ui/norm.ui"
#else
#define UI_FILE "src/melt.ui"
#endif

#include "element.hh"
#include "model.hh"
#include "phase.hh"
#include "update.hh"
#endif

#include "gui.hh"

//-----------------------------------------------------------------------------//
void on_save_as()
{
	Gtk::FileChooserDialog dialog("Please choose a file",
	                              Gtk::FILE_CHOOSER_ACTION_SAVE);
	dialog.set_transient_for(*main_win);

	//Add response buttons the the dialog:
	dialog.add_button("_Cancel", Gtk::RESPONSE_CANCEL);
	dialog.add_button("_Save", Gtk::RESPONSE_OK);

	//Add filters, so that only certain file types can be selected:

	Glib::RefPtr<Gtk::FileFilter> filter_text = Gtk::FileFilter::create();
	filter_text->set_name("Text files");
	filter_text->add_mime_type("text/plain");
	dialog.add_filter(filter_text);

	Glib::RefPtr<Gtk::FileFilter> filter_any = Gtk::FileFilter::create();
	filter_any->set_name("Any files");
	filter_any->add_pattern("*");
	dialog.add_filter(filter_any);

	//Show the dialog and wait for a user response:
	int result = dialog.run();

	//Handle the response:
	if (result == Gtk::RESPONSE_OK)
	{
		file = dialog.get_filename();
		do_save();
	}
}
