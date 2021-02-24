/*
 * melt.cc
 * Copyright (C) 2019 Kent G. Budge <kgb@kgbudge.com>
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

#include "Melt_Model.hh"


//-----------------------------------------------------------------------------//
Phase State::melt() const
{
	using namespace std;
    // Initial guess is that all fusible phases are fully melted. 
	// We will then see what should crystallize out.

	double XP[P_END];
	fill(XP, XP+P_END, 0.0);
	double Z[E_END];
	fill(Z, Z+E_END, 0.0);

	for (unsigned i=0; i<E_END; ++i)
	{
		double const x = X_[i];
		unsigned const p = ph_[i];
		if (x>0.0)
	    {
			Phase const &phase = phase_[p];
			unsigned const N = phase.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned const z = phase.z[j];
				unsigned const n = phase.n[j];
				Z[z] += x*n;
			}
			if (!is_fusible_[p])
			{
				XP[p] = x;
			}
        }
	}
	
	Melt_Model model(Z);

	return model.minimize_Gf(XP);
}

