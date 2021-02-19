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
                          unsigned const cphase[P_END],
                          double p[M_END],
                          double e) const
{
	using namespace std;
	
	double x[NP_];
	copy(X, X+NP_, x);
	for (unsigned i=0; i<n; ++i)
	{
		unsigned j = cphase[i];
		x[j] += e*p[i];
	}
	return Gf(x);
}
