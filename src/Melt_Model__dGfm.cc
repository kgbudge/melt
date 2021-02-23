/*
 * Melt_Model__dGfm.cc
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
#include <numeric>

//-----------------------------------------------------------------------------//
/*! Compute the derivative of melt free energy with the amount of a solid phase
 * extracted from the melt.
 * 
 * \param XP Number of moles of each solid phase.
 * 
 * \param xm Current melt composition, previously computed for XP.
 *
 * \param Gf0 Total free energy for the specified melt state, previously computed
 * for XP.
 * 
 * \param p Solid phase for which the derivative of total free energy is desired.
 * 
 * \return The derivate
 */
double Melt_Model::dGfm(double const cXP[P_END], 
                       double const xm[E_END],
                       double const Gf0, 
                       unsigned const p) const
{
	using namespace std;
	
	if (is_fusible_[p])
	{
		// See how much of the phase of interest can be crystallized from melt
		Phase const &ph = phase_[p];  
		unsigned const N = ph.nz;
		double nmax = std::numeric_limits<double>::max();
		for (int j=0; j<N; ++j)
		{
			unsigned z = ph.z[j];
			nmax = min(nmax, xm[z]/ph.n[j]);
		}
		double Gfp, Gfm, hh;
		if (nmax>1.0e-9) // We can crystallize a significant quantity.
		{
			double h = min(nmax, 0.01);
			double xmp[E_END];
			copy(xm_, xm_+E_END, xmp);		
			for (int j=0; j<N; ++j)
			{
				unsigned z = ph.z[j];
				xmp[z] -= h*ph.n[j];
			}
			Gfp = this->Gfm(xmp);
			hh = h;
		}
		else
		{
			hh = 0;
			Gfp = Gf0;
		}
		// We can always melt -- for purposes of calculating derivative, we "borrow" solid phase.
		{
			double h = min(nmax, 0.01);
			double xmp[E_END];
			copy(xm_, xm_+E_END, xmp);		
			for (int j=0; j<N; ++j)
			{
				unsigned z = ph.z[j];
				xmp[z] += h*ph.n[j];
			}
			Gfm = this->Gfm(xmp);
			hh += h;
		}
		return (Gfp-Gfm)/hh;
	}
	else
	{
		return 0.0; // not fusible
	}
}
