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


#include "kgb/Assert.h"
#include "kgb/tostring.h"

#include "Melt_Model.hh"
#include "Model.hh"
#include "State.hh"
#include "gui.hh"
#include "Phase.hh"
#include "update.hh"

//-----------------------------------------------------------------------------//
void update()
{
	using namespace std;

	string name;
	
	double nSiO2;
	double nTiO2;
	double nAl2O3;
	double nFe2O3;
	double nFeO;
	double nMnO;
	double nMgO;
	double nCaO;
	double nNa2O;
	double nK2O;
	double nP2O5;
	double nS;
	double nCr2O3;
	double nZrO2;
	double nH2O;
	double nCO2;
	double nCl;

	double total = get_molar_composition(name,
	                      nSiO2,
	                      nTiO2,
	                      nAl2O3,
	                      nFe2O3,
	                      nFeO,
	                      nMnO,
	                      nMgO,
	                      nCaO,
	                      nNa2O,
	                      nK2O,
	                      nP2O5,
	                      nS,
	                      nCr2O3,
	                      nZrO2,
	                      nH2O,
	                      nCO2,
	                      nCl);

	if (total==0.0) return;

	// Get the starting temperature and pressure
	string text = entry_T->get_text();
	double T = atof(text.c_str()) + 273.15;
	text = entry_P->get_text();
	double P = atof(text.c_str());

	// Initialize the State global data.
	State::initialize_globals(T, P, 
		          	          nH2O, nCO2, nNa2O, nMgO, nAl2O3, nSiO2, nP2O5, nS,
	                          nCl, nK2O, nCaO, nTiO2, nCr2O3, nMnO, nFeO,
	                          nFe2O3, nZrO2);

    // Initial guess is that all fusible phases are fully melted. 
	// We will then see what should crystallize out to get an optimum melt (possibly empty)

	// Create a State with the initial composition. 
	State state("1", 
	            nH2O, nCO2, nNa2O, nMgO, nAl2O3, nSiO2, nP2O5, nS, nCl,
	            nK2O, nCaO, nTiO2, nCr2O3, nMnO, nFeO, nFe2O3, nZrO2);

	// Calculate the equilibrium state at T and P.
	state.update();

	double px = 0.0;
	double const *const state_X = state.X();
	int const *const state_ph = state.ph();
	auto const &state_phase = state.phase();
	for (unsigned i=0; i<E_END; ++i)
	{
		if (state_X[i]>px && state_ph[i]+1==State::NP())
		{
			px = state_X[i];
		}
	}
	if (px > 0.0)
	{
		Phase const &melt = state_phase[State::NP()-1];

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
		
		text_nS->set_text(tostring(0.0));
		text_nCr2O3->set_text(tostring(0.0));
		text_nZrO2->set_text(tostring(0.0));
		text_nCO2->set_text(tostring(0.0));
		text_nCl->set_text(tostring(0.0));
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


	CIPW_Columns m_Columns;
	Glib::RefPtr<Gtk::ListStore> list_store_CIPW = Gtk::ListStore::create(m_Columns);
	list_store_CIPW->set_sort_column(1, Gtk::SORT_DESCENDING);

	Gtk::TreeModel::iterator iter;
	Gtk::TreeModel::Row row;

	// entry_pO2->set_text(tostring(pO2));

	// Convert to volume fraction

	double const *const state_V = state.V();
	double V[E_END];
	double Vtot = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		V[i] = state_V[i];
		Vtot += V[i];
	}

	double rnorm = 100.0/Vtot;
	for (unsigned i=0; i<E_END; ++i)
	{
		if (V[i]>0.0)
		{
			unsigned const pi = state_ph[i];
			iter = list_store_CIPW->append();
			row = *iter;
			row[m_Columns.m_col_text] = state_phase[pi].name;
			row[m_Columns.m_col_number] = V[i]*rnorm;
		}
	}

	tree_view_CIPW->remove_all_columns();
	tree_view_CIPW->set_model(list_store_CIPW);
	tree_view_CIPW->append_column("Mineral", m_Columns.m_col_text);
	tree_view_CIPW->append_column("Vol %", m_Columns.m_col_number);

	IUGS_classify(state);
}
