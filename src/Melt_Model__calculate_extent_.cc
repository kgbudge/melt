/*
 * Melt_Model__calculate_extent_.cc
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

#include <limits>

//-----------------------------------------------------------------------------//
/*! Calculate how far a reaction can proceed.
 *
 * \param reaction Reaction whose extent should be calculated.
 * \param XP Current solid phase composition
 */ 
double Melt_Model::calculate_extent_(Reaction const &reaction, 
                                     double const XP[],
                                     double const xm[],
                                     double const dGf) const
{
	 using namespace std;

	if (dGf<0)
	{
		// reaction goes forwards
	    double Result = numeric_limits<double>::max();
		int nz = reaction.nz;
		for (int i=0; i<nz; ++i)
		{
			if (reaction.n[i]<0.0)
			{
				Result = min(Result, -XP[reaction.p[i]]/reaction.n[i]);
			}
		}
		Phase const &phase = phase_[reaction.i];
		nz = phase.nz;
		for (int i=0; i<nz; ++i)
		{
			Result = min(Result, xm[phase.z[i]]/phase.n[i]);
		}
		return Result;
	}
	else
	{
		// reaction goes backwards
	    double Result = numeric_limits<double>::max();
		int nz = reaction.nz;
		for (int i=0; i<nz; ++i)
		{
			if (reaction.n[i]>0.0)
			{
				Result = min(Result, XP[reaction.p[i]]/reaction.n[i]);
			}
		}
		return Result;
	}
}
