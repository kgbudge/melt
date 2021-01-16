/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * on_open.cc
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

#include <fstream>
#include <iostream>

#include "gui.hh"

//-----------------------------------------------------------------------------//
void on_open()
{
	using namespace std;
	
	Gtk::FileChooserDialog dialog("Please choose a file",
	                              Gtk::FILE_CHOOSER_ACTION_OPEN);
	dialog.set_transient_for(*main_win);

	//Add response buttons the the dialog:
	dialog.add_button("_Cancel", Gtk::RESPONSE_CANCEL);
	dialog.add_button("_Open", Gtk::RESPONSE_OK);

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
		ifstream in(file.c_str());

		// Read the file contents

		on_new();
		string token;
		in >> token;
		while (in)
		{
			if (token=="name")
			{
				name.clear();
				char c = in.get();
				while (in && isspace(c) && c != '\n')
				{
					c = in.get();
				}
				while (in && c != '\n')
				{
					name += c;
					c = in.get();
				}
				entry_name->set_text(name);
			}
			else if (token=="SiO2")
			{
				in >> token;
				entry_SiO2->set_text(token);
			}
			else if (token=="TiO2")
			{
				in >> token;
				entry_TiO2->set_text(token);
			}
			else if (token=="Al2O3")
			{
				in >> token;
				entry_Al2O3->set_text(token);
			}
			else if (token=="Fe2O3")
			{
				in >> token;
				entry_Fe2O3->set_text(token);
			}
			else if (token=="FeO")
			{
				in >> token;
				entry_FeO->set_text(token);
			}
			else if (token=="MnO")
			{
				in >> token;
				entry_MnO->set_text(token);
			}
			else if (token=="MgO")
			{
				in >> token;
				entry_MgO->set_text(token);
			}
			else if (token=="CaO")
			{
				in >> token;
				entry_CaO->set_text(token);
			}
			else if (token=="Na2O")
			{
				in >> token;
				entry_Na2O->set_text(token);
			}
			else if (token=="K2O")
			{
				in >> token;
				entry_K2O->set_text(token);
			}
			else if (token=="P2O5")
			{
				in >> token;
				entry_P2O5->set_text(token);
			}
			else if (token=="S")
			{
				in >> token;
				entry_S->set_text(token);
			}
			else if (token=="Cr2O3")
			{
				in >> token;
				entry_Cr2O3->set_text(token);
			}
			else if (token=="ZrO2")
			{
				in >> token;
				entry_ZrO2->set_text(token);
			}
			else if (token=="H2O")
			{
				in >> token;
				entry_H2O->set_text(token);
			}
			else if (token=="CO2")
			{
				in >> token;
				entry_CO2->set_text(token);
			}
			else if (token=="Cl")
			{
				in >> token;
				entry_Cl->set_text(token);
			}
			else
			{
				cerr << token << '?' << endl;
				char c = in.get();
				while (in && c != '\n')
				{
					c = in.get();
				}
			}
			in >> token;
		}

		button_byweight->set_active(true);
		button_oxygen_by_composition->set_active(true);
		update();
	}
}
