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

#include "update.hh"

#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"

#include "ds++/Assert.hh"

#include "melt.hh"
#include "Model.hh"

//-----------------------------------------------------------------------------//
using namespace std;

unsigned const N = E_END;


//-----------------------------------------------------------------------------//
void calculate_Gf(double T /* K */, 
                  double P /* kbar */, 
                  bool is_element_active[E_END],
                  vector<Phase> &phase,
                  vector<double> &Gf /* kJ */,
                  vector<double> &amu /* g */)
{
	Require(phase.size()>=P_END);

	unsigned const N = phase.size();
	Gf.resize(N);
	amu.resize(N);
//	cout << "Potentially active phases:" << endl;
	for (unsigned i=0; i<N; ++i)
	{
		if (i!= phase[i].index)
		{
			cout << "ERROR: phase index out of synch for " << phase[i].name << endl;
			exit(1);
		}
		Gf[i] = phase[i].model->Gf(phase[i], T, P);
		double amui = 0.0;
		unsigned const nz = phase[i].nz;
		for (unsigned j=0; j<nz; ++j)
		{
			amui += atomic_weight[phase[i].z[j]]*phase[i].n[j];
			if (!is_element_active[phase[i].z[j]])
			{
				phase[i].nz = 0; // turn off this phase
			}
		}
		amu[i] = amui;
		if (phase[i].nz>0)
		{
//			cout << phase[i].name << " " << fixed << setprecision(3) << (1000*Gf[i]/amu[i]) << " kj/kg" << endl;
		}
	}

	// Water model selection
 	if (P>0.2224 || T>674.096 || T>649.634*pow(P,0.0811546)*(1+P*(1.39936 + P*(-6.98999+P*14.9787))))
	{
		  Gf[P_H2O_LIQUID] = 1e5; // turn liquid model off
	}
	else
	{
		  Gf[P_WATER_VAPOR] = 1e5; // turn vapor model off
    }
}
	
//-----------------------------------------------------------------------------//
bool 
do_ladder_update(double const T, 
                 double const P, 
                 vector<Phase> const &phase,
                 vector<double> Gf,
                 bool const oxygen_specified, 
                 bool const oxygen_FMQ,
                 double &pO2,
                 State &state)
{
	unsigned const P_END = phase.size();
	double GfO2;
	if (oxygen_specified)
	{
		if (oxygen_FMQ)
		{
			GfO2 = 2*Gf[P_MAGNETITE]+3*Gf[P_QUARTZ]-3*Gf[P_FAYALITE];
			pO2 = phase[P_O2].model->P(phase[P_O2], T, GfO2);
		}
		else
		{
			GfO2 = phase[P_O2].model->Gf(phase[P_O2], T, pO2);
		}
	}
	
	// Allocate storage for linear algebra operations

	unsigned const N = E_END;
	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_matrix *V = gsl_matrix_alloc(N, N);
	gsl_vector *S = gsl_vector_alloc(N);
	gsl_vector *work = gsl_vector_alloc(N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);	
	gsl_vector *aGf = gsl_vector_alloc(P_END);
	gsl_permutation *permute = gsl_permutation_alloc(P_END);
	double Ainv[N][N];
	double left[N];

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
		
		gsl_matrix_set_zero(A);
		gsl_vector_set_zero(b);

		if (oxygen_specified)
		{
			Gf[P_O2] = GfO2;
		}

		cout << "Active phases, state " << state.name <<  ':' << endl;
		for (unsigned i=0; i<N; ++i)
		{
			unsigned const pi = state.p[i];
			
			if (state.is_element_active[i]) 
			{
				cout << "  " << phase[pi].name << "  " << state.x[i] 
				    << " mol " << Gf[pi] << " kJ/mol" << endl;
			}
			
			gsl_vector_set(b, i, -Gf[pi]);
			unsigned const nz = phase[pi].nz;
			if (nz>0)
			{
			for (unsigned j=0; j<nz; j++)
			{
				gsl_matrix_set(A, i, phase[pi].z[j], phase[pi].n[j]);
			}
			}
			else
			{
				  gsl_matrix_set(A, i, i, 1.0);
			}
		}

		gsl_linalg_SV_decomp (A, V, S, work);		
		gsl_linalg_SV_solve (A, V, S, b, x);

		// The previous calculation should never be singular. We should get a
		// back activity for each element. This is the set of element activities
		// that would make the activity of every active phase zero -- a backwards
		// version of the usual free energy calculation where the element activities
		// are zero and this gives us a nonzero (and presumably negative!) free 
		// energy for each active phase. When we then use these nonzero element
		// activities to calculate new free energies for the library of phases,
		// any phase with a negative free energy is one that can be produced
		// by some reaction between existing phases to lower the energy further.

		cout << "Element activities for these phases:" << endl;
		for (unsigned i=0; i<E_END; ++i)
		{
			state.element_activity[i] = gsl_vector_get(x, i);
			if (state.is_element_active[i])
			{
				cout << "  " << element_name[i] << " " << state.element_activity[i] << " kJ/mol" << endl;
			}
		}

		// Look for candidate phases with negative backwards free energies. These
		// are phases that could be produced by reactions between existing phases
		// such that free energy is further lowered. However, we may find that
		// the reaction has nowhere to go -- the incoming phases are depleted.
		// This amounts to a ladder failure.
		cout << "New phase candidate activities:" << endl;

		bool found_candidate = false;
		for (unsigned i=0; i<P_END; ++i)
		{
			if (phase[i].nz>0)
			{
				double aGfi = Gf[i];
				double mol = 0.0;
				for (unsigned j=0; j<phase[i].nz; j++)
				{
					aGfi += phase[i].n[j]*state.element_activity[phase[i].z[j]];
					mol += phase[i].n[j];
				}
				gsl_vector_set(aGf, i, aGfi/mol);
				if (aGfi<-1.0e-9)
				{
					cout << "  " << phase[i].name << "  " << aGfi << " kJ/mol " << endl;
					found_candidate = true;
				}
			}
			else
			{
				gsl_vector_set(aGf, i, 10000.);
			}
		}
		if (!found_candidate)
		{
			success = true;
			goto DONE;
		}

		// Sort the free energy per mole atoms to find the best candidates for a new
		// phase. 
		gsl_sort_vector_index(permute, aGf);

		// Find a composition basis for the current composition. This is used
		// to construct the unique reaction between existing phases that
		// produces a candidate phase.

		for (unsigned i=0; i<N; i++) 
		{
			for (unsigned j=0; j<N; j++)
			{
				double sum = 0.0;
				for (unsigned k=0; k<N; ++k)
				{
					double sk = gsl_vector_get(S, k);
					if (fabs(sk)>1.0e-10)
					{
						sum += gsl_matrix_get(V, i, k)*gsl_matrix_get(A, j, k)/sk;
					}
				}
				Ainv[i][j] = sum;
			}
		}

		// Look through candidates for a new phase that lowers free energy.
		bool found = false;
		for (unsigned ii=0; ii<P_END; ++ii)
		{
			unsigned const ip = gsl_permutation_get(permute, ii);
			if (phase[ip].nz == 0 || oxygen_specified && ip == P_O2)
				continue;

			double aGfi = gsl_vector_get(aGf, ip); 
			if (fabs(aGfi)<1.0e-9)
			{
				break;
			}
			cout << "Considering phase " << phase[ip].name << ':' << endl;

			// Construct the reaction that produces the phase.

			fill(left, left+N, 0.0);
			for (unsigned i=0; i<phase[ip].nz; ++i)
			{
				unsigned z = phase[ip].z[i];
				double n = phase[ip].n[i];
				for (unsigned j=0; j<N; ++j)
				{
					left[j] += n*Ainv[z][j];
				}
			}

			// Calculate the reaction parameter, r, indicating how much of the 
			// reaction can take place.

			double r = numeric_limits<double>::max();
			double rGf = Gf[ip];
			for (unsigned i=0; i<N; ++i)
			{
				if (left[i]>1.0e-10)
				{
					r = min(r, state.x[i]/left[i]);
				}
				rGf -= left[i]*Gf[state.p[i]];
			}

			// We carry out the reaction even if the reaction parameter is zero.
			// Conceptually, we're replacing an infinitesimal amount of an old
			// phase with a new lower energy phase -- though both show zero
			// moles in the sample. This allows us to work our way down to
			// a real reaction that minimizes energy.
			cout << "Performing reaction ";
			found = true;
			bool first = true;
			for (unsigned i=0; i<N; ++i)
			{
				if (left[i]>1.0e-10)
				{
					if (!first) 
					{
						cout << " + ";
					}
					first = false;
					if (fabs(left[i]-1.0)>1.0e-10)
					{
						cout << setprecision(4) << left[i];
					}
					cout << phase[state.p[i]].name;
				}
			}
			cout << " -> " << phase[ip].name;
			
			for (unsigned i=0; i<N; ++i)
			{
				if (left[i]<-1.0e-10)
				{
					cout << " + ";
					first = false;
					if (fabs(left[i]+1.0)>1.0e-10)
					{
						cout << setprecision(4) << -left[i];
					}
					cout << phase[state.p[i]].name;
				}
			}
			cout << endl;

			// Now carry out the reaction, noting which reagents are exhausted.

			unsigned n = 0;
 	        unsigned p[E_END]; 
			for (unsigned i=0; i<N; ++i)
			{
				state.x[i] -= r*left[i];
				if (fabs(state.x[i])<1.0e-10)
				{
					if (left[i]>1.0e-9)
					{
						cout << "Reaction depletes " << phase[state.p[i]].name << endl;
						p[n++] = i;
					}
					state.x[i] = 0.0;
				}
			}
			// n should not be zero
			if (n==1)
			{
				state.x[p[0]] = r;	
				state.p[p[0]] = ip;
				break;
			}
			else
			{
				// we've exhausted two or more phases simultaneously. Split the state.
				vector<State> states(n, state);
				double Gfn = 0.0;
				for (unsigned i=0; i<E_END; ++i)
				{
					Gfn += state.x[i]*Gf[state.p[i]];
				}
				cout << "  Starting subladder Gf = " << Gfn << " kJ" << endl;
				cout << "Substate ladder:" << endl;
				for (unsigned i=0; i<n; ++i)
				{
					states[i].x[p[i]] = r;	
					states[i].p[p[i]] = ip;
					states[i].name += '.' + to_string(i+1);
					cout << "Substate " << states[i].name << ':' << endl;
					
					bool result = do_ladder_update(T, 
					                             P, 
					                             phase,
					                             Gf,
					                             oxygen_specified, 
					                             oxygen_FMQ,
					                             pO2,
					                             states[i]);

					cout << "Comparing branch " << states[i].name << " against best so far " << endl;
					double Gfnn = 0.0;
					for (unsigned j=0; j<E_END; ++j)
					{
						Gfnn += states[i].x[j]*Gf[states[i].p[j]];
					}
					cout << "  Gf = " << Gfnn << " kJ" << endl;
					if (Gfnn<Gfn)
					{
						cout << "  State " << states[i].name << " is new best." << endl;	
						Gfn = Gfnn;
						state = states[i];
					}
					else
					{
						cout << "Best unchanged." << endl;
					}
				}
				goto DONE;
			}
		}
		if (!found)
		{
			cout << "No reaction found." << endl;
			goto DONE;
		}
	}


	DONE:
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(aGf);
	gsl_permutation_free(permute);

	if (true || !oxygen_specified)
	{
		pO2 = phase[P_O2].model->P(phase[P_O2], T, -2*state.element_activity[E_O]);
	}
		
	// Convert to volume fraction

	double Vtot = 0.0;
	for (unsigned i=0; i<N; ++i)
	{
		unsigned const pi = state.p[i];
        state.V[i] = state.x[i]*phase[pi].model->volume(phase[pi], T, P);
		Vtot += state.V[i];
	}
		
	double rnorm = 100.0/Vtot;
	for (unsigned i=0; i<N; ++i)
	{
		state.V[i] *= rnorm;
	}

	return success;
}

//-----------------------------------------------------------------------------//
void update_state(double const T, 
                  double const P, 
                  vector<Phase> &phase,
                  bool const oxygen_specified, 
                  bool const oxygen_FMQ,
                  double &pO2,
                  State &state)
{
	// Determine which elements are active 

	for (unsigned e=0; e<E_END; ++e)
	{
		state.is_element_active[e] = (state.x[e]>0.0);
	}
	cout << "Active elemental phases:" << endl;
	for (unsigned e=0; e<E_END; ++e)
	{
		if (state.is_element_active[e])
		{
			cout << " " << phase[state.p[e]].name << endl;
		}
	}
	
	vector<double> Gf, amu;
	// Calculate the free energy of formation per mole and molecular weight for all possible phases
	calculate_Gf(T, P, state.is_element_active, phase, Gf, amu);

	// Calculate the free energy of the sample
	double Gftot = 0.0;
	double Mtot = 0.0;
	for (unsigned e=0; e<E_END; ++e)
	{
		Gftot += state.x[e]*Gf[state.p[e]];
		Mtot += state.x[e]*amu[state.p[e]];
	}
	cout << "Initial free energy of formation = " << (1000*Gftot/Mtot) << " kJ" << endl;

	for (;;)
	{
		cout << endl << "Ladder search for minimum free energy:" << endl;
		auto status = do_ladder_update(T, 
		                               P, 
		                               phase, 
		                               Gf, 
		                               oxygen_specified, 
		                               oxygen_FMQ, 
		                               pO2, 
		                               state);

		cout << "Active phases:" << endl;
		for (unsigned e=0; e<E_END; ++e)
		{
			if (state.x[e]>0.0)
			{
				cout << " " << phase[state.p[e]].name << " = " << state.x[e] << endl;
			}
		}

		Gftot = 0.0;
		Mtot = 0.0;
		for (unsigned e=0; e<E_END; ++e)
		{
			Gftot += state.x[e]*Gf[state.p[e]];
			Mtot += state.x[e]*amu[state.p[e]];
		}
		cout << "Free energy of formation of sample = " << fixed << setprecision(3) << Gftot << " kJ" << endl;

		for (unsigned i=0; i<5; ++i) // No iteration for now -- will implement multiple melt later
		{
			Phase new_phase;
			cout << "  Starting free energy for this melt step: " << Gftot << " kJ" << endl;
			double Geu  = melt(T, P, phase, Gf, state, new_phase);
			
			if (Geu < Gftot - 1e-9)
			{
				Gf.push_back(new_phase.Hf0);
				phase.push_back(new_phase);

				do_ladder_update(T, 
				                 P, 
				                 phase, 
				                 Gf, 
				                 oxygen_specified, 
				                 oxygen_FMQ, 
				                 pO2, 
				                 state);

				Gftot = 0.0;
				for (unsigned e=0; e<E_END; ++e)
				{
					Gftot += state.x[e]*Gf[state.p[e]];
				}
				cout << "    Free energy of formation = " << fixed << setprecision(3) << (1000*Gftot/Mtot) << " kJ/kg" << endl;
			}
			// else we've done all the melting we can
			else
			{
				return;
			}
		}
		return;
	}
}

