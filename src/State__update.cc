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

#include "State.hh"

#include <iomanip>
#include <iostream>
#include <numeric>

#include "Model.hh"

//-----------------------------------------------------------------------------//
using namespace std;

//-----------------------------------------------------------------------------//
void State::update()
{
	using namespace std;

	// Calculate the free energy of the sample
	double Gftot = 0.0;
	for (unsigned e=0; e<E_END; ++e)
	{
		Gftot += X_[e]*Gf_[ph_[e]];
	}
	cout << "Initial free energy of formation = " << Gftot << " kJ" << endl;

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

	cout << "  Starting free energy for this melt step: " << Gftot << " kJ" << endl;
	Phase new_phase = melt();

	Gf_[NP_-1] = new_phase.Hf0;
	phase_[NP_-1] = new_phase;

	do_ladder_update();

	Gftot = 0.0;
	for (unsigned e=0; e<E_END; ++e)
	{
		Gftot += X_[e]*Gf_[ph_[e]];
	}
	cout << "    Free energy of formation = " << fixed 
		<< setprecision(3) << Gftot << " kJ" << endl;
	
	for (unsigned i=0; i<E_END; ++i)
    {
		Phase const &phase = phase_[ph_[i]];
		V_[i] = X_[i]*phase.model->volume(phase, T_, P_);
	}
	return;
}

