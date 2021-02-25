/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * State__do_ladder_update.cc
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

#include "State.hh"

#include <iomanip>
#include <iostream>
#include <limits>

#include "gsl/gsl_sort_vector.h"

#include "Model.hh"

//-----------------------------------------------------------------------------//
bool 
State::do_ladder_update()
{
	using namespace std;
	
	// Allocate storage for linear algebra operations

	double Ainv[E_END][E_END];
	double left[E_END];

	bool success;
	
	// now iterate to find composition

	for(;;)
	{
		// Find absolute free energy of elements corresponding to current
		// composition. This is the set of free energy values for the elements
		// that makes the free energy of every active phase zero.

		// There are always as many active phases as active elements. How does this
		// work if there is a degenerate composition, in which a phase completely
		// consumes two elements? This is the reason for the ladder. If a reaction
		// simultaneously depletes more than one existing phase to create the
		// new lower energy phase, we fork our calculation, with a different
		// depleted phase replace by the new phase in each fork. We pick the
		// fork that minimizes total free energy.

		// Build a linear system and solve for element free energies
		
		gsl_matrix_set_zero(gsl_A_);
		gsl_vector_set_zero(gsl_b_);

		for (unsigned i=0; i<E_END; ++i)
		{
			unsigned const pi = ph_[i];			
			gsl_vector_set(gsl_b_, i, -Gf_[pi]);
			Phase const &phase = phase_[pi];
			unsigned const nz = phase.nz;
			for (unsigned j=0; j<nz; j++)
			{
				gsl_matrix_set(gsl_A_, i, phase.z[j], phase.n[j]);
			}
		}

		gsl_linalg_SV_decomp (gsl_A_, gsl_V_, gsl_S_, gsl_work_);		
		gsl_linalg_SV_solve (gsl_A_, gsl_V_, gsl_S_, gsl_b_, gsl_x_);

		// The previous calculation should never be singular. We should get a
		// back activity for each element. This is the set of element activities
		// that would make the activity of every active phase zero -- a backwards
		// version of the usual free energy calculation where the element activities
		// are zero and this gives us a nonzero (and presumably negative!) free 
		// energy for each active phase. When we then use these nonzero element
		// activities to calculate new free energies for the library of phases,
		// any phase with a negative free energy is one that can be produced
		// by some reaction between existing phases to lower the energy further.

		for (unsigned i=0; i<E_END; ++i)
		{
			element_activity_[i] = gsl_vector_get(gsl_x_, i);
		}

		// Look for candidate phases with negative backwards free energies. These
		// are phases that could be produced by reactions between existing phases
		// such that free energy is further lowered. However, we may find that
		// the reaction has nowhere to go -- the incoming phases are depleted.
		// This amounts to a ladder failure.
//		cout << "New phase candidate activities:" << endl;

		bool found_candidate = false;
		for (unsigned i=0; i<NP_; ++i)
		{
			double aGfi = Gf_[i];
			double mol = 0.0;
			for (unsigned j=0; j<phase_[i].nz; j++)
			{
				aGfi += phase_[i].n[j]*element_activity_[phase_[i].z[j]];
				mol += phase_[i].n[j];
			}
			gsl_vector_set(gsl_aGf_, i, aGfi/mol);
			if (aGfi<-1.0e-9)
			{
//				cout << "  " << phase_[i].name << "  " << aGfi << " kJ/mol " << endl;
				found_candidate = true;
			}
		}
		if (!found_candidate)
		{
			success = true;
			goto DONE;
		}

		// Sort the free energy per mole atoms to find the best candidates for a new
		// phase. 
		gsl_sort_vector_index(gsl_permute_, gsl_aGf_);

		// Find a composition basis for the current composition. This is used
		// to construct the unique reaction between existing phases that
		// produces a candidate phase.

		for (unsigned i=0; i<E_END; i++) 
		{
			for (unsigned j=0; j<E_END; j++)
			{
				double sum = 0.0;
				for (unsigned k=0; k<E_END; ++k)
				{
					double sk = gsl_vector_get(gsl_S_, k);
					if (fabs(sk)>1.0e-10)
					{
						sum += gsl_matrix_get(gsl_V_, i, k)*gsl_matrix_get(gsl_A_, j, k)/sk;
					}
				}
				Ainv[i][j] = sum;
			}
		}

		// Look through candidates for a new phase that lowers free energy.
		bool found = false;
		for (unsigned ii=0; ii<NP_; ++ii)
		{
			unsigned const ip = gsl_permutation_get(gsl_permute_, ii);
			if (phase_[ip].nz == 0)
				continue;

			double aGfi = gsl_vector_get(gsl_aGf_, ip); 
			if (fabs(aGfi)<1.0e-9)
			{
				break;
			}
//			cout << "Considering phase " << phase_[ip].name << ':' << endl;

			// Construct the reaction that produces the phase.

			fill(left, left+E_END, 0.0);
			for (unsigned i=0; i<phase_[ip].nz; ++i)
			{
				unsigned z = phase_[ip].z[i];
				double n = phase_[ip].n[i];
				for (unsigned j=0; j<E_END; ++j)
				{
					left[j] += n*Ainv[z][j];
				}
			}

			// Calculate the reaction parameter, r, indicating how much of the 
			// reaction can take place.

			double r = numeric_limits<double>::max();
			double rGf = Gf_[ip];
			for (unsigned i=0; i<E_END; ++i)
			{
				if (left[i]>1.0e-10)
				{
					r = min(r, X_[i]/left[i]);
				}
				rGf -= left[i]*Gf_[ph_[i]];
			}

			// We carry out the reaction even if the reaction parameter is zero.
			// Conceptually, we're replacing an infinitesimal amount of an old
			// phase with a new lower energy phase -- though both show zero
			// moles in the sample. This allows us to work our way down to
			// a real reaction that minimizes energy.
			//cout << "Performing reaction ";
			found = true;
			bool first = true;
			for (unsigned i=0; i<E_END; ++i)
			{
				if (left[i]>1.0e-10)
				{
					if (!first) 
					{
						//cout << " + ";
					}
					first = false;
					if (fabs(left[i]-1.0)>1.0e-10)
					{
						//cout << setprecision(4) << left[i];
					}
					//cout << phase_[ph_[i]].name;
				}
			}
			//cout << " -> " << phase_[ip].name;
			
			for (unsigned i=0; i<E_END; ++i)
			{
				if (left[i]<-1.0e-10)
				{
					//cout << " + ";
					first = false;
					if (fabs(left[i]+1.0)>1.0e-10)
					{
						//cout << setprecision(4) << -left[i];
					}
					//cout << phase_[ph_[i]].name;
				}
			}
			//cout << endl;

			// Now carry out the reaction, noting which reagents are exhausted.

			unsigned n = 0;
 	        unsigned pd[E_END]; 
			for (unsigned i=0; i<E_END; ++i)
			{
				X_[i] -= r*left[i];
				if (fabs(X_[i])<1.0e-10)
				{
					if (left[i]>1.0e-9)
					{
						//cout << "Reaction depletes " << phase_[ph_[i]].name << endl;
						pd[n++] = i;
					}
					X_[i] = 0.0;
				}
			}
			// n should not be zero
			// By allowing null reactions, we eliminate need for ladder update!
			X_[pd[0]] = r;	
			ph_[pd[0]] = ip;
			break;
		}
		if (!found)
		{
			goto DONE;
		}
	}


	DONE:
		
	// Convert to volume fraction

	double Vtot = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		unsigned const pi = ph_[i];
		if (pi+1 != NP_)
		{
          V_[i] = X_[i]*phase_[pi].model->volume(phase_[pi], T_, P_);
		  Vtot += V_[i];
		}
		else
		{
			V_[i] = X_[i]*phase_[pi].V;
			Vtot += V_[i];
		}
	}
		
	double rnorm = 100.0/Vtot;
	for (unsigned i=0; i<E_END; ++i)
	{
		V_[i] *= rnorm;
	}

	return success;
}
