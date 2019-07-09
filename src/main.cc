/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cc
 * Copyright (C) 2015 Kent G. Budge <kgb@kgbudge.com>
 * 
 * norm is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * norm is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
/* #define UI_FILE PACKAGE_DATA_DIR"/ui/norm.ui" */
#define UI_FILE "src/melt.ui"

#include "phase.hh"

//-----------------------------------------------------------------------------//
class CIPW_Columns : public Gtk::TreeModelColumnRecord
{
public:

  CIPW_Columns()
    { add(m_col_text); add(m_col_number); }

  Gtk::TreeModelColumn<Glib::ustring> m_col_text;
  Gtk::TreeModelColumn<double> m_col_number;
};

//-----------------------------------------------------------------------------//
using namespace std;

Gtk::Window* main_win;

Gtk::Entry *entry_name;

Gtk::Entry *entry_SiO2;
Gtk::Entry *entry_TiO2;
Gtk::Entry *entry_Al2O3;
Gtk::Entry *entry_Fe2O3;
Gtk::Entry *entry_FeO;
Gtk::Entry *entry_MnO;
Gtk::Entry *entry_MgO;
Gtk::Entry *entry_CaO;
Gtk::Entry *entry_Na2O;
Gtk::Entry *entry_K2O;
Gtk::Entry *entry_P2O5;
Gtk::Entry *entry_S;
Gtk::Entry *entry_Cr2O3;
Gtk::Entry *entry_ZrO2;
Gtk::Entry *entry_H2O;
Gtk::Entry *entry_CO2;
Gtk::Entry *entry_Cl;
Gtk::Entry *entry_T;
Gtk::Entry *entry_P;


Gtk::Label *text_nSiO2;
Gtk::Label *text_nTiO2;
Gtk::Label *text_nAl2O3;
Gtk::Label *text_nFe2O3;
Gtk::Label *text_nFeO;
Gtk::Label *text_nMnO;
Gtk::Label *text_nMgO;
Gtk::Label *text_nCaO;
Gtk::Label *text_nNa2O;
Gtk::Label *text_nK2O;
Gtk::Label *text_nP2O5;
Gtk::Label *text_nS;
Gtk::Label *text_nCr2O3;
Gtk::Label *text_nZrO2;
Gtk::Label *text_nH2O;
Gtk::Label *text_nCO2;
Gtk::Label *text_nCl;

Gtk::Label *text_total;

Gtk::Label *text_ibc_family;
Gtk::Label *text_tas;
Gtk::Label *text_iugs;

Gtk::RadioButton *button_byweight;
Gtk::RadioButton *button_oxygen_by_composition;
Gtk::RadioButton *button_oxygen_specified;
Gtk::RadioButton *button_oxygen_FMQ;
Gtk::Entry *entry_pO2;

Gtk::TreeView *tree_view_CIPW;

double SiO2;
double TiO2;
double Al2O3;
double Fe2O3;
double FeO;
double MnO;
double MgO;
double CaO;
double Na2O;
double K2O;
double P2O5;
double S;
double Cr2O3;
double ZrO2;
double H2O;
double CO2;
double Cl;

string file;
string name;

void update(double const T, 
            double const P, 
            vector<Phase> phase,
            bool const oxygen_specified, 
            bool const oxygen_FMQ,
            double &pO2,
	        double Np[E_END],
	        unsigned p[E_END],
            double volume[E_END]);

//-----------------------------------------------------------------------------//
void update()
{
	name = entry_name->get_text();

	string text = entry_SiO2->get_text();
	SiO2 = atof(text.c_str());
	text = entry_TiO2->get_text();
	TiO2 = atof(text.c_str());
	text = entry_Al2O3->get_text();
	Al2O3 = atof(text.c_str());
	text = entry_Fe2O3->get_text();
	Fe2O3 = atof(text.c_str());
	text = entry_FeO->get_text();
	FeO = atof(text.c_str());
	text = entry_MnO->get_text();
	MnO = atof(text.c_str());
	text = entry_MgO->get_text();
	MgO = atof(text.c_str());
	text = entry_CaO->get_text();
	CaO = atof(text.c_str());
	text = entry_Na2O->get_text();
	Na2O = atof(text.c_str());
	text = entry_K2O->get_text();
	K2O = atof(text.c_str());
	text = entry_P2O5->get_text();
	P2O5 = atof(text.c_str());
	text = entry_S->get_text();
	S = atof(text.c_str());
	text = entry_Cr2O3->get_text();
	Cr2O3 = atof(text.c_str());
	text = entry_ZrO2->get_text();
	ZrO2 = atof(text.c_str());
	text = entry_H2O->get_text();
	H2O = atof(text.c_str());
	text = entry_CO2->get_text();
	CO2 = atof(text.c_str());
	text = entry_Cl->get_text();
	Cl = atof(text.c_str());
	text = entry_T->get_text();
	double T = atof(text.c_str());
	text = entry_P->get_text();
	double P = atof(text.c_str());

	double total = SiO2 + TiO2 + Al2O3 + Fe2O3 + FeO + MnO + MgO + CaO + Na2O + K2O + P2O5 + S + Cr2O3 + ZrO2 + H2O + CO2 + Cl;
	text_total->set_text(tostring(total));
	
	double nSiO2 = SiO2;
	double nTiO2 = TiO2;
	double nAl2O3 = Al2O3;
	double nFe2O3 = Fe2O3;
	double nFeO = FeO;
	double nMnO = MnO;
	double nMgO = MgO;
	double nCaO = CaO;
	double nNa2O = Na2O;
	double nK2O = K2O;
	double nP2O5 = P2O5;
	double nS = S;
	double nCr2O3 = Cr2O3;
	double nZrO2 = ZrO2;
	double nH2O = H2O;
	double nCO2 = CO2;
	double nCl = Cl;
	
	// Calculate molar composition if selected

	if (button_byweight->get_active())
	{
		nSiO2 /= 60.085;
		nTiO2 /= 79.899;
		nAl2O3 /= 101.961;
		nFe2O3 /= 159.692;
		nFeO /= 71.85;
		nMnO /= 70.94;
		nMgO /= 40.311;
		nCaO /= 56.08;
		nNa2O /= 61.98;
		nK2O /= 94.20;
		nP2O5 /= 109.95;
		nS /= 32.06;
		nCr2O3 /= 223.84;
		nZrO2 /= 123.22;
		nH2O /= 18.02;
		nCO2 /= 44.01;
		nCl /= 35.45;
	}
	
	total = nSiO2 + nTiO2 + nAl2O3 + nFe2O3 + nFeO + nMnO + nMgO + nCaO + nNa2O + nK2O + nP2O5 + nS + nCr2O3 + nZrO2 + nH2O + nCO2 + nCl;

	if (total==0.0) return;

	double rnorm = 100.0/total;
	nSiO2 *= rnorm;
	text_nSiO2->set_text(tostring(nSiO2));
	nTiO2 *= rnorm;
	text_nTiO2->set_text(tostring(nTiO2));
	nAl2O3 *= rnorm;
	text_nAl2O3->set_text(tostring(nAl2O3));
	nFe2O3 *= rnorm;
	text_nFe2O3->set_text(tostring(nFe2O3));
	nFeO *= rnorm;
	text_nFeO->set_text(tostring(nFeO));
	nMnO *= rnorm;
	text_nMnO->set_text(tostring(nMnO));
	nMgO *= rnorm;
	text_nMgO->set_text(tostring(nMgO));
	nCaO *= rnorm;
	text_nCaO->set_text(tostring(nCaO));
	nNa2O *= rnorm;
	text_nNa2O->set_text(tostring(nNa2O));
	nK2O *= rnorm;
	text_nK2O->set_text(tostring(nK2O));
	nP2O5 *= rnorm;
	text_nP2O5->set_text(tostring(nP2O5));
	nS *= rnorm;
	text_nS->set_text(tostring(nS));
	nCr2O3 *= rnorm;
	text_nCr2O3->set_text(tostring(nCr2O3));
	nZrO2 *= rnorm;
	text_nZrO2->set_text(tostring(nZrO2));
	nH2O *= rnorm;
	text_nH2O->set_text(tostring(nH2O));
	nCO2 *= rnorm;
	text_nCO2->set_text(tostring(nCO2));
	nCl *= rnorm;
	text_nCl->set_text(tostring(nCl));

	// Initial composition
	double Np[E_END];
	unsigned p[E_END];
	double volume[E_END];

	double Gf[P_END];
	initialize_Gf(T, P, Gf);

	bool oxygen_specified = !button_oxygen_by_composition->get_active();
	bool oxygen_FMQ;
	double pO2, GfO2;
	if (oxygen_specified)
	{
		if (button_oxygen_specified->get_active())
		{
			text = entry_pO2->get_text();
			pO2 = atof(text.c_str());
			GfO2 = phase[P_O2].model->Gf(phase[P_O2], T, pO2);
			oxygen_FMQ = false;
		}
		else
		{
			Check(button_oxygen_FMQ->get_active());
			GfO2 = 2*Gf[P_MAGNETITE]+3*Gf[P_SiO2_QUARTZ]-3*Gf[P_FAYALITE];
			pO2 = phase[P_O2].model->P(phase[P_O2], T, GfO2);
			oxygen_FMQ = true;
		}
	}

	p[E_H] = P_H2O_LIQUID;
	Np[E_H] = nH2O;

	p[E_C] = P_CO2;
	Np[E_C] = nCO2;

	if (!oxygen_specified)
	{
	  p[E_O] = P_Fe2O3;
  	  Np[E_O] = nFe2O3;	
	}
	else
	{
		p[E_O] = P_O2;
		Np[E_O] = numeric_limits<double>::max();
	}

	p[E_NA] = P_Na2O;
	Np[E_NA] = nNa2O - 0.5*nCl;

	p[E_MG] = P_MgO;
	Np[E_MG] = nMgO;

	p[E_AL] = P_Al2O3;
	Np[E_AL] = nAl2O3;

	p[E_SI] = P_SiO2_QUARTZ;
	Np[E_SI] = nSiO2;

	p[E_S] = P_S;
	Np[E_S] = nS;

	p[E_CL] = P_HALITE;
	Np[E_CL] = nCl;

	p[E_K] = P_K2O;
	Np[E_K] = nK2O;

	p[E_CA] = P_CaO;
	Np[E_CA] = nCaO;

	p[E_TI] = P_TiO2;
	Np[E_TI] = nTiO2;

	p[E_CR] = P_Cr2O3;
	Np[E_CR] = nCr2O3;

	p[E_MN] = P_MnO;
	Np[E_MN] = nMnO;

	p[E_FE] = P_FeO;
	Np[E_FE] = nFeO;
	if (oxygen_specified)
	{
		Np[E_FE] += 2*nFe2O3;
	}

	p[E_ZR] = P_ZrO2;
	Np[E_ZR] = nZrO2;
	
    update(T, 
           P, 
           vector<Phase>(phase, phase+P_END),
           oxygen_specified, 
           oxygen_FMQ,
           pO2,
	       Np,
	       p,
           volume);

	CIPW_Columns m_Columns;
	Glib::RefPtr<Gtk::ListStore> list_store_CIPW = Gtk::ListStore::create(m_Columns);
	list_store_CIPW->set_sort_column(1, Gtk::SORT_DESCENDING);

	Gtk::TreeModel::iterator iter;
	Gtk::TreeModel::Row row;

	entry_pO2->set_text(tostring(pO2));
		
	// Convert to volume fraction

	double Vtot = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		unsigned const pi = p[i];
        volume[i] = Np[i]*phase[pi].model->volume(phase[pi], T, P);
		Vtot += volume[i];
	}
		
	rnorm = 100.0/Vtot;
	for (unsigned i=0; i<E_END; ++i)
	{
		volume[i] *= rnorm;
		if (volume[i]>0.0)
		{
			unsigned const pi = p[i];
			iter = list_store_CIPW->append();
			row = *iter;
			row[m_Columns.m_col_text] = phase[pi].name;
			row[m_Columns.m_col_number] = volume[i];
		}
	}

	tree_view_CIPW->remove_all_columns();
	tree_view_CIPW->set_model(list_store_CIPW);
	tree_view_CIPW->append_column("Mineral", m_Columns.m_col_text);
	tree_view_CIPW->append_column("Vol %", m_Columns.m_col_number);
}

//-----------------------------------------------------------------------------//
void on_new()
{
	file = "";
	name = "";
	entry_name->set_text("");
	entry_SiO2->set_text("0.00");
	entry_TiO2->set_text("0.00");
	entry_Al2O3->set_text("0.00");
	entry_Fe2O3->set_text("0.00");
	entry_FeO->set_text("0.00");
	entry_MnO->set_text("0.00");
	entry_MgO->set_text("0.00");
	entry_CaO->set_text("0.00");
	entry_Na2O->set_text("0.00");
	entry_K2O->set_text("0.00");
	entry_P2O5->set_text("0.00");
	entry_S->set_text("0.00");
	entry_Cr2O3->set_text("0.00");
	entry_ZrO2->set_text("0.00");
	entry_H2O->set_text("0.00");
	entry_CO2->set_text("0.00");
	entry_Cl->set_text("0.00");
	update();
}

//-----------------------------------------------------------------------------//
void on_open()
{
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

//-----------------------------------------------------------------------------//
void do_save()
{
	ofstream out(file.c_str());
	if (!out)
	{
		Gtk::MessageDialog dialog(*main_win, "Save Fail");
		dialog.set_secondary_text("Could not write save file " + file);

		dialog.run();
	}
	else
	{
		out << "name " << name << endl;
		out << "SiO2 " << SiO2 << endl;
		out << "TiO2 " << TiO2 << endl;
		out << "Al2O3 " << Al2O3 << endl;
		out << "Fe2O3 " << Fe2O3 << endl;
		out << "FeO " << FeO << endl;
		out << "MnO " << MnO << endl;
		out << "MgO " << MgO << endl;
		out << "CaO " << CaO << endl;
		out << "Na2O " << Na2O << endl;
		out << "K2O " << K2O << endl;
		out << "P2O5 " << P2O5 << endl;
		out << "S " << S << endl;
		out << "Cr2O3 " << Cr2O3 << endl;
		out << "ZrO2 " << ZrO2 << endl;
		out << "H2O " << H2O << endl;
		out << "CO2 " << CO2 << endl;
		out << "Cl " << Cl << endl;
	}
}

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

//-----------------------------------------------------------------------------//
void on_save()
{
	if (file=="")
	{
		on_save_as();
	}
	else
	{
		do_save();
	}
}

//-----------------------------------------------------------------------------//
void on_quit()
{
	main_win->hide();
}

//-----------------------------------------------------------------------------//
void normalize()
{
	double total = SiO2 + TiO2 + Al2O3 + Fe2O3 + FeO + MnO + MgO + CaO + Na2O + K2O + P2O5 + S + Cr2O3 + ZrO2 + H2O + CO2 + Cl;
	if (total>0.0)
	{
		double rnorm = 100.0/total;
		entry_SiO2->set_text(tostring(SiO2*rnorm));
		entry_TiO2->set_text(tostring(TiO2*rnorm));
		entry_Al2O3->set_text(tostring(Al2O3*rnorm));
		entry_Fe2O3->set_text(tostring(Fe2O3*rnorm));
		entry_FeO->set_text(tostring(FeO*rnorm));
		entry_MnO->set_text(tostring(MnO*rnorm));
		entry_MgO->set_text(tostring(MgO*rnorm));
		entry_CaO->set_text(tostring(CaO*rnorm));
		entry_Na2O->set_text(tostring(Na2O*rnorm));
		entry_K2O->set_text(tostring(K2O*rnorm));
		entry_P2O5->set_text(tostring(P2O5*rnorm));
		entry_S->set_text(tostring(S*rnorm));
		entry_Cr2O3->set_text(tostring(Cr2O3*rnorm));
		entry_ZrO2->set_text(tostring(ZrO2*rnorm));
		entry_H2O->set_text(tostring(H2O*rnorm));
		entry_CO2->set_text(tostring(CO2*rnorm));
		entry_Cl->set_text(tostring(Cl*rnorm));
		update();
	}
}

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
		entry_name->signal_activate().connect(sigc::ptr_fun(update));

		builder->get_widget("entry_SiO2", entry_SiO2);
		entry_SiO2->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_TiO2", entry_TiO2);
		entry_TiO2->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_Al2O3", entry_Al2O3);
		entry_Al2O3->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_Fe2O3", entry_Fe2O3);
		entry_Fe2O3->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_FeO", entry_FeO);
		entry_FeO->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_MnO", entry_MnO);
		entry_MnO->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_MgO", entry_MgO);
		entry_MgO->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_CaO", entry_CaO);
		entry_CaO->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_Na2O", entry_Na2O);
		entry_Na2O->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_K2O", entry_K2O);
		entry_K2O->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_P2O5", entry_P2O5);
		entry_P2O5->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_S", entry_S);
		entry_S->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_Cr2O3", entry_Cr2O3);
		entry_Cr2O3->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_ZrO2", entry_ZrO2);
		entry_ZrO2->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_H2O", entry_H2O);
		entry_H2O->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_CO2", entry_CO2);
		entry_CO2->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_Cl", entry_Cl);
		entry_Cl->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_T", entry_T);
		entry_T->signal_activate().connect(sigc::ptr_fun(update));
		builder->get_widget("entry_P", entry_P);
		entry_P->signal_activate().connect(sigc::ptr_fun(update));

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
		button->signal_clicked().connect(sigc::ptr_fun(normalize));
		builder->get_widget("button_recalculate", button);
		button->signal_clicked().connect(sigc::ptr_fun(update));
		builder->get_widget("button_bymol", button);
		button->signal_activate().connect(sigc::ptr_fun(update));

		builder->get_widget("button_byweight", button_byweight);
		button_byweight->signal_activate().connect(sigc::ptr_fun(update));

		builder->get_widget("button_oxygen_by_composition", button_oxygen_by_composition);
		builder->get_widget("button_oxygen_specified", button_oxygen_specified);
		builder->get_widget("button_oxygen_FMQ", button_oxygen_FMQ);

		builder->get_widget("entry_pO2", entry_pO2);
		entry_pO2->signal_activate().connect(sigc::ptr_fun(update));

		builder->get_widget("text_total", text_total);

		builder->get_widget("text_ibc_family", text_ibc_family);

		builder->get_widget("text_tas", text_tas);

		builder->get_widget("text_iugs", text_iugs);

		builder->get_widget("tree_view_CIPW", tree_view_CIPW);

		// Go for it
		kit.run(*main_win);
	}
	return 0;
}


