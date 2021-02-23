/*
 * Melt_Model__Gfmelt.cc
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
double Melt_Model::Gfmelt(double const X[P_END],
                          unsigned const n,
                          Reaction const cphase[P_END],
                          double p[M_END],
                          double e)
{
	using namespace std;
	
	double x[NP_];
	copy(X, X+NP_, x);
	double Result = 0.0;
	for (unsigned i=0; i<n; ++i)
	{
		auto const &r = cphase[i];
		unsigned const N = r.nz;
		for (unsigned j=0; j<N; ++j)
		{		
			unsigned const ph = r.p[j];
		    x[ph] += e*p[i]*r.n[j];
	  	    Result += x[ph]*Gf_[ph];
		}
	}
	double xm[E_END];
	compute_current_melt_composition_(x, xm);
	Result += Gfm(xm);
	return Result;
}
