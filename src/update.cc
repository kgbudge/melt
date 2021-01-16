/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * update.cc
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
Gtk::Entry *entry_over_rho;
Gtk::Entry *entry_depth;


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
Gtk::RadioButton *button_molar_melt;
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
#endif

#include "kgb/Assert.h"
#include "kgb/tostring.h"

#include "Model.hh"
#include "State.hh"
#include "composition.hh"
#include "gui.hh"
#include "phase.hh"
#include "update.hh"

//-----------------------------------------------------------------------------//
void update()
{
	using namespace std;
	
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
	double T = atof(text.c_str()) + 273.15;
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
	State state;
	double Gf[P_END];

	bool oxygen_specified = !button_oxygen_by_composition->get_active();
	bool oxygen_FMQ;
	double pO2, GfO2;
	if (oxygen_specified)
	{
		if (button_oxygen_specified->get_active())
		{
			text = entry_pO2->get_text();
			pO2 = atof(text.c_str());
			oxygen_FMQ = false;
		}
		else
		{
			Check(button_oxygen_FMQ->get_active());
			oxygen_FMQ = true;
		}
	}

	state.p[E_H] = P_H2;
	state.x[E_H] = nH2O;
	double nO2 = 0.5*nH2O;

	state.p[E_C] = P_GRAPHITE;
	state.x[E_C] = nCO2;
	nO2 += nCO2;

	state.p[E_NA] = P_Na;
	state.x[E_NA] = 2*nNa2O;
	nO2 += 0.5*nNa2O;

	state.p[E_MG] = P_Mg;
	state.x[E_MG] = nMgO;
	nO2 += 0.5*nMgO;

	state.p[E_AL] = P_Al;
	state.x[E_AL] = 2*nAl2O3;
	nO2 += 1.5*nAl2O3;

	state.p[E_SI] = P_Si;
	state.x[E_SI] = nSiO2;
	nO2 += nSiO2;

	state.p[E_P] = P_P4;
	state.x[E_P] = 0.5*nP2O5;
	nO2 += 2.5*nP2O5;

	state.p[E_S] = P_S;
	state.x[E_S] = nS;

	state.p[E_CL] = P_Cl2;
	state.x[E_CL] = 0.5*nCl;

	state.p[E_K] = P_K;
	state.x[E_K] = 2*nK2O;
	nO2 += 0.5*nK2O;

	state.p[E_CA] = P_Ca;
	state.x[E_CA] = nCaO;
	nO2 += 0.5*nCaO;

	state.p[E_TI] = P_Ti;
	state.x[E_TI] = nTiO2;
	nO2 += nTiO2;

	state.p[E_CR] = P_Cr;
	state.x[E_CR] = 2*nCr2O3;
	nO2 += 1.5*nCr2O3;

	state.p[E_MN] = P_Mn;
	state.x[E_MN] = nMnO;
	nO2 += 0.5*nMnO;

	state.p[E_FE] = P_Fe;
	state.x[E_FE] = nFeO + 2*nFe2O3;
	nO2 += 0.5*nFeO + 1.5*nFe2O3;

	state.p[E_ZR] = P_Zr;
	state.x[E_ZR] = nZrO2;
	nO2 += nZrO2;

	state.p[E_O] = P_O2;
	state.x[E_O] = nO2;

	vector<Phase> phase(::phase, ::phase+P_END);

	state.name="1";

	update_state(T, 
	             P, 
	             phase,
	             oxygen_specified, 
	             oxygen_FMQ,
	             pO2,
	             state);
	
	if (button_molar_melt->get_active())
	{
		int pm = -1;
		double px = 0.0;
		for (unsigned i=0; i<E_END; ++i)
		{
			if (state.x[i]>px && state.p[i]>=P_END)
			{
				pm = state.p[i];
				px = state.x[i];
			}
		}
		if (pm != -1)
		{
			Phase const &melt = phase[pm];

			unsigned const N = melt.nz;
			double mH2O = 0;
			double mSiO2 = 0;
			double mAl2O3 = 0;
			double mMgO = 0;
			double mFeO = 0;
			double mCaO = 0;
			double mNa2O = 0;
			double mK2O = 0;
			for (unsigned i=0; i<N; ++i)
			{
				double const ni = melt.n[i];
				switch (melt.z[i])
				{
					case E_H:
						mH2O += 0.5*ni;
						break;

					case E_NA:
						mNa2O += 0.5*ni;
						break;

					case E_MG:
						mMgO += ni;
						break;

					case E_AL:
						mAl2O3 += 0.5*ni;
						break;

					case E_SI:
						mSiO2 += ni;
						break;

					case E_K:
						mK2O += 0.5*ni;
						break;

					case E_CA:
						mCaO += ni;
						break;

					case E_FE:
						mFeO +=  ni;
						break;

					default:
						break;
				}
			}

			double tot = mH2O + mSiO2 + mAl2O3 + mMgO + mFeO + mCaO + mK2O + mNa2O;
			double rtot = 100/tot;

			text_nH2O->set_text(tostring(mH2O*rtot));
			text_nSiO2->set_text(tostring(mSiO2*rtot));
			text_nAl2O3->set_text(tostring(mAl2O3*rtot));
			text_nMgO->set_text(tostring(mMgO*rtot));
			text_nFeO->set_text(tostring(mFeO*rtot));
			text_nCaO->set_text(tostring(mCaO*rtot));
			text_nK2O->set_text(tostring(mK2O*rtot));
			text_nNa2O->set_text(tostring(mNa2O*rtot));

			text_nTiO2->set_text(tostring(0.0));
			text_nFe2O3->set_text(tostring(0.0));
			text_nMnO->set_text(tostring(0.0));
			text_nP2O5->set_text(tostring(0.0));
		}
		else
		{
			text_nH2O->set_text("--");
			text_nSiO2->set_text("--");
			text_nAl2O3->set_text("--");
			text_nMgO->set_text("--");
			text_nFeO->set_text("--");
			text_nCaO->set_text("--");
			text_nK2O->set_text("--");
			text_nNa2O->set_text("--");

			text_nTiO2->set_text("--");
			text_nFe2O3->set_text("--");
			text_nMnO->set_text("--");
			text_nP2O5->set_text("--");
		}

	}
	
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
		if (state.x[i]>0.0)
		{
			unsigned const pi = state.p[i];
			if (pi<P_END)
			{
				state.V[i] = state.x[i]*phase[pi].model->volume(phase[pi], T, P);
			}
			else
			{
				state.V[i] = phase[pi].V;
			}
			Vtot += state.V[i];
		}
		else
		{
			state.V[i] = 0.0;
		}
	}
	//	cout << Vtot/100 << endl;
	
	rnorm = 100.0/Vtot;
	for (unsigned i=0; i<E_END; ++i)
	{
		state.V[i] *= rnorm;
		if (state.V[i]>0.0)
		{
			unsigned const pi = state.p[i];
			iter = list_store_CIPW->append();
			row = *iter;
			row[m_Columns.m_col_text] = phase[pi].name;
			row[m_Columns.m_col_number] = state.V[i];
		}
	}

	tree_view_CIPW->remove_all_columns();
	tree_view_CIPW->set_model(list_store_CIPW);
	tree_view_CIPW->append_column("Mineral", m_Columns.m_col_text);
	tree_view_CIPW->append_column("Vol %", m_Columns.m_col_number);
}
