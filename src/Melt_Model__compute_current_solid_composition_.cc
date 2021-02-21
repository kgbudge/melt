/*
 * Melt_Model__compute_current_solid_composition.cc
 * Copyright (C) 2021 Kent G. Budge <kgb@kgbudge.com>
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

#include "Melt_Model.hh"

//-----------------------------------------------------------------------------//
void Melt_Model::compute_current_solid_composition_(double const XP[])
{
	using namespace std;
	
	fill(xs_, xs_+E_END, 0.0);
	for (unsigned i=0; i<NP_; ++i)
	{
		double const x = XP[i];
		if (x>0.0)
		{
			Phase const &phase = phase_[i];  
			unsigned const N = phase.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned z = phase.z[j];
				xs_[z] += x*phase.n[j];
			}
		}
	}
}
