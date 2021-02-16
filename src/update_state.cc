/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * update.cc
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

#include "update.hh"

#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>
#include <cmath>

#include "ds++/Assert.hh"

#include "Melt_Model.hh"
#include "Model.hh"

//-----------------------------------------------------------------------------//
using namespace std;

unsigned const N = E_END;


//-----------------------------------------------------------------------------//
void calculate_Gf(double T /* K */, 
                  double P /* kbar */, 
                  bool is_element_active[E_END],
                  vector<Phase> &phase,
                  vector<double> &Gf /* kJ */,
                  vector<double> &amu /* g */)
{
	Require(phase.size()>=P_END);

	unsigned const N = phase.size();
	Gf.resize(N);
	amu.resize(N);
//	cout << "Potentially active phases:" << endl;
	for (unsigned i=0; i<N; ++i)
	{
		if (i!= phase[i].index)
		{
			cout << "ERROR: phase index out of synch for " << phase[i].name << endl;
			exit(1);
		}
		Gf[i] = phase[i].model->Gf(phase[i], T, P);
		double amui = 0.0;
		unsigned const nz = phase[i].nz;
		for (unsigned j=0; j<nz; ++j)
		{
			amui += atomic_weight[phase[i].z[j]]*phase[i].n[j];
			if (!is_element_active[phase[i].z[j]])
			{
				phase[i].nz = 0; // turn off this phase
			}
		}
		amu[i] = amui;
		if (phase[i].nz>0)
		{
//			cout << phase[i].name << " " << fixed << setprecision(3) << (1000*Gf[i]/amu[i]) << " kj/kg" << endl;
		}
	}

	// Water model selection
 	if (P>0.2224 || T>674.096 || T>649.634*pow(P,0.0811546)*(1+P*(1.39936 + P*(-6.98999+P*14.9787))))
	{
		  Gf[P_H2O_LIQUID] = 1e5; // turn liquid model off
	}
	else
	{
		  Gf[P_WATER_VAPOR] = 1e5; // turn vapor model off
    }
}

//-----------------------------------------------------------------------------//
void State::update()
{
	unsigned const NPL = phase_.size();

	// Calculate the free energy of the sample
	double Gftot = 0.0;
	for (unsigned e=0; e<E_END; ++e)
	{
		Gftot += X_[e]*Gf_[ph_[e]];
	}
	cout << "Initial free energy of formation = " << Gftot << " kJ" << endl;

	for (;;)
	{
		cout << endl << "Ladder search for minimum free energy:" << endl;
		auto status = do_ladder_update();

		cout << "Active phases:" << endl;
		for (unsigned e=0; e<E_END; ++e)
		{
			if (X_[e]>0.0)
			{
				cout << " " << phase_[ph_[e]].name << " = " << X_[e] << endl;
			}
		}

		Gftot = 0.0;
		for (unsigned e=0; e<E_END; ++e)
		{
			Gftot += X_[e]*Gf_[ph_[e]];
		}
		cout << "Free energy of formation of sample = " << fixed << setprecision(3) << Gftot << " kJ" << endl;

		for (unsigned i=0; i<5; ++i) // No iteration for now -- will implement multiple melt later
		{
			cout << "  Starting free energy for this melt step: " << Gftot << " kJ" << endl;
			Phase new_phase = melt();
			
			if (new_phase.nz>0)
			{
				Gf_.push_back(new_phase.Hf0);
				phase_.push_back(new_phase);

				do_ladder_update();

				Gftot = 0.0;
				for (unsigned e=0; e<E_END; ++e)
				{
					Gftot += X_[e]*Gf_[ph_[e]];
				}
				cout << "    Free energy of formation = " << fixed 
					<< setprecision(3) << Gftot << " kJ" << endl;
			}
			// else we've done all the melting we can
			else
			{
				for (unsigned i=0; i<E_END; ++i)
				{
					Phase const &phase = phase_[ph_[i]];
					V_[i] = X_[i]*phase.model->volume(phase, T_, P_);
				}
				return;
			}
		}
		return;
	}
}

