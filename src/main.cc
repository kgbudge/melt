/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cc
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

#include <iostream>

/* For testing propose use the local (not installed) ui file */
#ifndef NDEBUG
#define UI_FILE PACKAGE_DATA_DIR"/ui/norm.ui"
#else
#define UI_FILE "src/melt.ui"
#endif

#define EXTERN
#include "composition.hh"
#include "gui.hh"

//-----------------------------------------------------------------------------//
int
main (int argc, char *argv[])
{
	Gtk::Main kit(argc, argv);


	//Load the Glade file and instiate its widgets:
	Glib::RefPtr<Gtk::Builder> builder;
	try
	{
		builder = Gtk::Builder::create_from_file(UI_FILE);
	}
	catch (const Glib::FileError & ex)
	{
		std::cerr << ex.what() << std::endl;
		return 1;
	}
	builder->get_widget("main_window", main_win);


	if (main_win)
	{
		// Hook up all the actions
		Gtk::ImageMenuItem *menu_item;
		builder->get_widget("menu_file_new", menu_item);
		menu_item->signal_activate().connect(sigc::ptr_fun(on_new)); 
		builder->get_widget("menu_file_open", menu_item);
		menu_item->signal_activate().connect(sigc::ptr_fun(on_open)); 
		builder->get_widget("menu_file_save", menu_item);
		menu_item->signal_activate().connect(sigc::ptr_fun(on_save)); 
		builder->get_widget("menu_file_save_as", menu_item);
		menu_item->signal_activate().connect(sigc::ptr_fun(on_save_as)); 
		builder->get_widget("menu_file_quit", menu_item);
		menu_item->signal_activate().connect(sigc::ptr_fun(on_quit)); 
				
		builder->get_widget("entry_name", entry_name);
		entry_name->signal_activate().connect(sigc::ptr_fun(on_activate));

		builder->get_widget("entry_SiO2", entry_SiO2);
		entry_SiO2->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_TiO2", entry_TiO2);
		entry_TiO2->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_Al2O3", entry_Al2O3);
		entry_Al2O3->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_Fe2O3", entry_Fe2O3);
		entry_Fe2O3->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_FeO", entry_FeO);
		entry_FeO->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_MnO", entry_MnO);
		entry_MnO->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_MgO", entry_MgO);
		entry_MgO->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_CaO", entry_CaO);
		entry_CaO->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_Na2O", entry_Na2O);
		entry_Na2O->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_K2O", entry_K2O);
		entry_K2O->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_P2O5", entry_P2O5);
		entry_P2O5->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_S", entry_S);
		entry_S->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_Cr2O3", entry_Cr2O3);
		entry_Cr2O3->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_ZrO2", entry_ZrO2);
		entry_ZrO2->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_H2O", entry_H2O);
		entry_H2O->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_CO2", entry_CO2);
		entry_CO2->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_Cl", entry_Cl);
		entry_Cl->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_T", entry_T);
		entry_T->signal_activate().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("entry_P", entry_P);
		entry_P->signal_activate().connect(sigc::ptr_fun(on_activate_P));
		builder->get_widget("entry_over_rho", entry_over_rho);
		builder->get_widget("entry_depth", entry_depth);
		entry_depth->signal_activate().connect(sigc::ptr_fun(on_activate_depth));

		builder->get_widget("text_nSiO2", text_nSiO2);
		builder->get_widget("text_nTiO2", text_nTiO2);
		builder->get_widget("text_nAl2O3", text_nAl2O3);
		builder->get_widget("text_nFe2O3", text_nFe2O3);
		builder->get_widget("text_nFeO", text_nFeO);
		builder->get_widget("text_nMnO", text_nMnO);
		builder->get_widget("text_nMgO", text_nMgO);
		builder->get_widget("text_nCaO", text_nCaO);
		builder->get_widget("text_nNa2O", text_nNa2O);
		builder->get_widget("text_nK2O", text_nK2O);
		builder->get_widget("text_nP2O5", text_nP2O5);
		builder->get_widget("text_nS", text_nS);
		builder->get_widget("text_nCr2O3", text_nCr2O3);
		builder->get_widget("text_nZrO2", text_nZrO2);
		builder->get_widget("text_nH2O", text_nH2O);
		builder->get_widget("text_nCO2", text_nCO2);
		builder->get_widget("text_nCl", text_nCl);

		Gtk::Button *button;
		builder->get_widget("button_norm", button);
		button->signal_clicked().connect(sigc::ptr_fun(on_normalize));
		builder->get_widget("button_recalculate", button);
		button->signal_clicked().connect(sigc::ptr_fun(on_activate));
		builder->get_widget("button_bymol", button);
		button->signal_activate().connect(sigc::ptr_fun(on_activate));

		builder->get_widget("button_byweight", button_byweight);
		button_byweight->signal_activate().connect(sigc::ptr_fun(on_activate));

		builder->get_widget("button_oxygen_by_composition", button_oxygen_by_composition);
		builder->get_widget("button_oxygen_specified", button_oxygen_specified);
		builder->get_widget("button_oxygen_FMQ", button_oxygen_FMQ);

		builder->get_widget("button_molar_melt", button_molar_melt);

		builder->get_widget("entry_pO2", entry_pO2);
		entry_pO2->signal_activate().connect(sigc::ptr_fun(on_activate));

		builder->get_widget("text_total", text_total);

		builder->get_widget("text_iugs", text_iugs);

		builder->get_widget("tree_view_CIPW", tree_view_CIPW);

		// Go for it
		kit.run(*main_win);
	}
	return 0;
}


