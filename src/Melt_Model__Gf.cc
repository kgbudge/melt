/*
 * Melt_Model__Gf.cc
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
double Melt_Model::Gf(double const XP[P_END]) const
{
	using namespace std;
	
	unsigned const NM = NP_;

	double xm[E_END];
	for (unsigned i=0; i<E_END; ++i)
	{
		xm[i] = Z_[i];
	}
	// Modify for amount of each fusible phase crystallized out
	double Result = Gfr_;
	for (unsigned i=0; i<NM; ++i)
	{
		double const x = XP[i];
		Check(x>=0.0);
		if (x>0.0)
		{
			Phase const &phase = phase_[i];  
			unsigned const N = phase.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned z = phase.z[j];
				xm[z] -= x*phase.n[j];
				xm[z] = max(0.0, xm[z]);
			}
			Result += x*Gf_[i];
		}
	}
    Result += this->Gfm(xm);

	return Result;
}
