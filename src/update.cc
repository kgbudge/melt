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

#include "Model.hh"
#include "State.hh"
#include "gui.hh"
#include "phase.hh"
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
			string text = entry_pO2->get_text();
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

	string text = entry_T->get_text();
	double T = atof(text.c_str()) + 273.15;
	text = entry_P->get_text();
	double P = atof(text.c_str());

	update_state(T, 
	             P, 
	             phase,
	             oxygen_specified, 
	             oxygen_FMQ,
	             pO2,
	             state);
	
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
	
	double rnorm = 100.0/Vtot;
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

	IUGS_classify(state);
}
