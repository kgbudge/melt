/*
 * Melt_Model__mimimize_Gf.cc
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

#include <iomanip>
#include <iostream>
#include <limits>
#include <set>

#include "Model.hh"

//-----------------------------------------------------------------------------//
/*! Find the mixture of solid phases and a melt phase that minimized total Gibbs
 * free energy.
 * 
 * \param[in,out] XP Contains initial guess of number of moles of each solid 
 * phase. This is compressed to include only the Melt::NP_ active phases 
 * specified by Melt::p_. On return, contains number of moles of each solid 
 * phase that minimizes free energy.
 *
 * \return The melt phase that minimizes free energy.
 */
Phase Melt_Model::minimize_Gf(double XP[P_END])
{
	using namespace std;

	// Restart loop. This discards the search space and begins building up a
	// new search space for the minimum of free energy.
	for (;;)
	{
		Reaction cphase[E_END];
		unsigned NMP = 0;
		update_current_state_(XP);
		
		// Main loop. Here we build up a search space for the minimum of free energy.
		for (;;)
		{
			// Calculate free energy of current estimate of sample state, and its
			// gradient with every solid phase in library.

			double Gfm = this->Gfm(xm_);
			double Gf = this->Gf(XP);
			cout << "Gfm = " << fixed << setprecision(3) << Gfm << endl;
			cout << "Gf = " << fixed << setprecision(3) << Gf << endl;
			double best = numeric_limits<double>::max();
			Reaction rbest;
			for (unsigned i=0; i<NP_; ++i)
			{
				if (is_fusible_[i])
				{
					Reaction r = construct_reaction_(XP, xm_, xs_, Gf, Gfm, i);
					if (r.extent*r.dGf0<best)
					{
						best = r.extent*r.dGf0;
						rbest = r;
					}
				}
			}

			// If nothing can change, we must be done.

			if (best>1.0e-9*cnorm_ || NMP+1==E_END)
			{
				break;
			}

			// Add a single phase that is our best guess of the next melt element

			cphase[NMP] = rbest;
			NMP++;

			// Minimize on this melt set.

			double revised_Gf = minimize_trial_set_(NMP, cphase, XP);
			update_current_state_(XP);

			double mtot = 0.0;
			for (unsigned i=0; i<E_END; ++i)
			{
				if (xm_[i]<1.0e-9*cnorm_)
				{
					xm_[i] = 0.0;
				}
				mtot += xm_[i];
			}	
			if (mtot<=0.0 || Gf - revised_Gf  < 1.0e-9*cnorm_)
			{
				break;
			}
		}
	}
		
	double V = 0.0;
	double xend[M_END];
	compute_melt_endmember_composition_(xm_, xend);
	for (unsigned i=0; i<M_END; ++i)
	{
		double xi = xend[i];
		if (i != E_O && xi>0.0)
		{
			Phase const &phase = ::phase[melt_endmember[i]];
			V += xi*phase.model->volume(phase, T_, P_);
		}
	}
		
	Phase Result;
	Result.index = 0;
	Result.name = "melt";
	Result.nz = 0;
	Result.V = V;
	for (unsigned i=0; i<E_END; ++i)
	{
		if (xm_[i]>1e-9)
		{
			Result.z[Result.nz] = i;
			Result.n[Result.nz] = xm_[i];
			Result.nz++;
		}
	}
	// If Result.nz is zero, this flags to client that there is no melt.

	Result.S0 = 0.0;
	Result.Hf0 = Gfm(xm_);
	Result.model = MELT;

	return Result;
}
