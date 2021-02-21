/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * State.cc
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

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"

#include "Assert.hh"

#include "Melt_Model.hh"
#include "Model.hh"
#include "Phase.hh"
#include "phase_enum.hh"

//-----------------------------------------------------------------------------//
State::State(std::string const &name, 
             double T,
             double P,
             double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
             double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
             double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
             double nFe2O3,  double nZrO2)

: name_(name), T_(T), P_(P)
{
	Require(T>0.0);
	Require(P>0.0);
	Require(nH2O>=0.0);
	Require(nCO2>=0.0);  
	Require(nNa2O>=0.0);
	Require(nMgO>=0.0);
	Require(nAl2O3>=0.0);
	Require(nSiO2>=0.0);
	Require(nP2O5>=0.0);  
	Require(nS>=0.0); 
	Require(nCl>=0.0); 
	Require(nK2O>=0.0);
	Require(nCaO>=0.0); 
	Require(nTiO2>=0.0);
	Require(nCr2O3>=0.0);  
	Require(nMnO>=0.0); 
	Require(nFeO>=0.0);
	Require(nFe2O3>=0.0); 
	Require(nZrO2>=0.0);
	
	using namespace std;

	ph_[E_H] = P_H2;
	X_[E_H] = nH2O;
	is_element_active_[E_H] = nH2O>0.0;
	double nO2 = 0.5*nH2O;

	ph_[E_C] = P_GRAPHITE;
	X_[E_C] = nCO2;
	is_element_active_[E_C] = nCO2>0.0;
	nO2 += nCO2;

	ph_[E_NA] = P_Na;
	X_[E_NA] = 2*nNa2O;
	is_element_active_[E_NA] = nNa2O>0.0;
	nO2 += 0.5*nNa2O;

	ph_[E_MG] = P_Mg;
	X_[E_MG] = nMgO;
	is_element_active_[E_MG] = nMgO>0.0;
	nO2 += 0.5*nMgO;

	ph_[E_AL] = P_Al;
	X_[E_AL] = 2*nAl2O3;
	is_element_active_[E_AL] = nAl2O3>0.0;
	nO2 += 1.5*nAl2O3;

	ph_[E_SI] = P_Si;
	X_[E_SI] = nSiO2;
	is_element_active_[E_SI] = nSiO2>0.0;
	nO2 += nSiO2;

	ph_[E_P] = P_P4;
	X_[E_P] = 0.5*nP2O5;
	is_element_active_[E_P] = nP2O5>0.0;
	nO2 += 2.5*nP2O5;

	ph_[E_S] = P_S;
	X_[E_S] = nS;
	is_element_active_[E_S] = nS>0.0;

	ph_[E_CL] = P_Cl2;
	X_[E_CL] = 0.5*nCl;
	is_element_active_[E_CL] = nCl>0.0;

	ph_[E_K] = P_K;
	X_[E_K] = 2*nK2O;
	is_element_active_[E_K] = nK2O;
	nO2 += 0.5*nK2O;

	ph_[E_CA] = P_Ca;
	X_[E_CA] = nCaO;
	nO2 += 0.5*nCaO;
	is_element_active_[E_CA] = nCaO>0.0;

	ph_[E_TI] = P_Ti;
	X_[E_TI] = nTiO2;
	nO2 += nTiO2;
	is_element_active_[E_TI] = nTiO2>0.0;

	ph_[E_CR] = P_Cr;
	X_[E_CR] = 2*nCr2O3;
	nO2 += 1.5*nCr2O3;
	is_element_active_[E_CR] = nCr2O3>0.0;

	ph_[E_MN] = P_Mn;
	X_[E_MN] = nMnO;
	nO2 += 0.5*nMnO;
	is_element_active_[E_MN] = nMnO>0.0;

	ph_[E_FE] = P_Fe;
	X_[E_FE] = nFeO + 2*nFe2O3;
	nO2 += 0.5*nFeO + 1.5*nFe2O3;
	is_element_active_[E_FE] = X_[E_FE]>0.0;

	ph_[E_ZR] = P_Zr;
	X_[E_ZR] = nZrO2;
	nO2 += nZrO2;
	is_element_active_[E_ZR] = nZrO2>0.0;

	ph_[E_O] = P_O2;
	X_[E_O] = nO2;
	is_element_active_[E_O] = nO2>0.0;

	compute_gf_();
}

//-----------------------------------------------------------------------------//
State::State(std::string const &name, 
             double T,
             double P,
             double const x0[E_END])

: name_(name), T_(T), P_(P)
{
	Require(T>0.0);
	Require(P>0.0);
	Require(std::count_if(x0, x0+E_END, [](double x){return x<0.0;}) == 0);
	
	using namespace std;

	copy(x0, x0+E_END, X_);
	for (unsigned i=0; i<E_END; ++i)
	{
		is_element_active_[i] = x0[i]>0.0;
	}

	ph_[E_H] = P_H2;
	X_[E_H] *= 0.5;
	ph_[E_C] = P_GRAPHITE;
	ph_[E_NA] = P_Na;
	ph_[E_MG] = P_Mg;
	ph_[E_AL] = P_Al;
	ph_[E_SI] = P_Si;
	ph_[E_P] = P_P4;
	X_[E_P] *= 0.25;
	ph_[E_S] = P_S;
	ph_[E_CL] = P_Cl2;
	X_[E_CL] *= 0.5;
	ph_[E_K] = P_K;
	ph_[E_CA] = P_Ca;
	ph_[E_TI] = P_Ti;
	ph_[E_CR] = P_Cr;
	ph_[E_MN] = P_Mn;
	ph_[E_FE] = P_Fe;
	ph_[E_ZR] = P_Zr;
	ph_[E_O] = P_O2;
	X_[E_O] *= 0.5;

	compute_gf_();
}
	
//-----------------------------------------------------------------------------//
bool 
State::do_ladder_update()
{
	using namespace std;
	
	unsigned const NPL = phase_.size();
	
	// Allocate storage for linear algebra operations

	gsl_matrix *A = gsl_matrix_alloc(E_END, E_END);
	gsl_matrix *g_V = gsl_matrix_alloc(E_END, E_END);
	gsl_vector *S = gsl_vector_alloc(E_END);
	gsl_vector *work = gsl_vector_alloc(E_END);
	gsl_vector *b = gsl_vector_alloc(E_END);
	gsl_vector *g_x = gsl_vector_alloc(E_END);	
	gsl_vector *aGf = gsl_vector_alloc(NPL);
	gsl_permutation *permute = gsl_permutation_alloc(NPL);
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
		
		gsl_matrix_set_zero(A);
		gsl_vector_set_zero(b);

		for (unsigned i=0; i<E_END; ++i)
		{
			unsigned const pi = ph_[i];			
			gsl_vector_set(b, i, -Gf_[pi]);
			Phase const &phase = phase_[pi];
			unsigned const nz = phase.nz;
			for (unsigned j=0; j<nz; j++)
			{
				gsl_matrix_set(A, i, phase.z[j], phase.n[j]);
			}
		}

		gsl_linalg_SV_decomp (A, g_V, S, work);		
		gsl_linalg_SV_solve (A, g_V, S, b, g_x);

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
			element_activity_[i] = gsl_vector_get(g_x, i);
		}

		// Look for candidate phases with negative backwards free energies. These
		// are phases that could be produced by reactions between existing phases
		// such that free energy is further lowered. However, we may find that
		// the reaction has nowhere to go -- the incoming phases are depleted.
		// This amounts to a ladder failure.
//		cout << "New phase candidate activities:" << endl;

		bool found_candidate = false;
		for (unsigned i=0; i<NPL; ++i)
		{
			double aGfi = Gf_[i];
			double mol = 0.0;
			for (unsigned j=0; j<phase_[i].nz; j++)
			{
				aGfi += phase_[i].n[j]*element_activity_[phase_[i].z[j]];
				mol += phase_[i].n[j];
			}
			gsl_vector_set(aGf, i, aGfi/mol);
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
		gsl_sort_vector_index(permute, aGf);

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
					double sk = gsl_vector_get(S, k);
					if (fabs(sk)>1.0e-10)
					{
						sum += gsl_matrix_get(g_V, i, k)*gsl_matrix_get(A, j, k)/sk;
					}
				}
				Ainv[i][j] = sum;
			}
		}

		// Look through candidates for a new phase that lowers free energy.
		bool found = false;
		for (unsigned ii=0; ii<NPL; ++ii)
		{
			unsigned const ip = gsl_permutation_get(permute, ii);
			if (phase_[ip].nz == 0)
				continue;

			double aGfi = gsl_vector_get(aGf, ip); 
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
			cout << "Performing reaction ";
			found = true;
			bool first = true;
			for (unsigned i=0; i<E_END; ++i)
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
					cout << phase_[ph_[i]].name;
				}
			}
			cout << " -> " << phase_[ip].name;
			
			for (unsigned i=0; i<E_END; ++i)
			{
				if (left[i]<-1.0e-10)
				{
					cout << " + ";
					first = false;
					if (fabs(left[i]+1.0)>1.0e-10)
					{
						cout << setprecision(4) << -left[i];
					}
					cout << phase_[ph_[i]].name;
				}
			}			cout << endl;

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
						cout << "Reaction depletes " << phase_[ph_[i]].name << endl;
						pd[n++] = i;
					}
					X_[i] = 0.0;
				}
			}
			// n should not be zero
			if (n==1)
			{
				X_[pd[0]] = r;	
				ph_[pd[0]] = ip;
				break;
			}
			else
			{
				// we've exhausted two or more phases simultaneously. Split the state.
				vector<State> states(n, *this);
				double Gfn = r*Gf_[ip];
				for (unsigned i=0; i<E_END; ++i)
				{
					Gfn += X_[i]*Gf_[ph_[i]];
				}
				for (unsigned i=0; i<n; ++i)
				{
					states[i].X_[pd[i]] = r;	
					states[i].ph_[pd[i]] = ip;
					states[i].name_ += '.' + to_string(i+1);
					bool result = states[i].do_ladder_update();
					double Gfnn = 0.0;
					for (unsigned j=0; j<E_END; ++j)
					{
						Gfnn += states[i].X_[j]*Gf_[states[i].ph_[j]];
					}
					if (Gfnn<Gfn)
					{
						Gfn = Gfnn;
						*this = states[i];
					}
				}
				goto DONE;
			}
		}
		if (!found)
		{
			goto DONE;
		}
	}


	DONE:
	gsl_matrix_free(A);
	gsl_matrix_free(g_V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_vector_free(b);
	gsl_vector_free(g_x);
	gsl_vector_free(aGf);
	gsl_permutation_free(permute);
		
	// Convert to volume fraction

	double Vtot = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		unsigned const pi = ph_[i];
        V_[i] = X_[i]*phase_[pi].model->volume(phase_[pi], T_, P_);
		Vtot += V_[i];
	}
		
	double rnorm = 100.0/Vtot;
	for (unsigned i=0; i<E_END; ++i)
	{
		V_[i] *= rnorm;
	}

	return success;
}

void State::compute_gf_()
{
	using namespace std;
	
	phase_.reserve(P_END);
	Gf_.reserve(P_END);

	// Basic phases are always active.
	for (unsigned i=0; i<E_END; ++i)
	{
		Phase const &phase = ::phase[i];
		if (i!= phase.index)
		{
			cout << "ERROR: phase index out of synch for " << phase.name << endl;
			exit(1);
		}
		phase_.push_back(phase);
		Gf_.push_back(phase.model->Gf(phase, T_, P_));
	}

	bool no_liquid_water = P_>0.2224 || T_>674.096 ||
		T_>649.634*pow(P_,0.0811546)*(1+P_*(1.39936 + P_*(-6.98999+P_*14.9787)));

	for (unsigned i=E_END; i<P_END; i++)
	{
		Phase const &phase = ::phase[i];
		if (i!= phase.index)
		{
			cout << "ERROR: phase index out of synch for " << phase.name << endl;
			exit(1);
		}
		unsigned const N = phase.nz;
		bool active = true;
		for (unsigned j=0; j<N; ++j)
		{
			if (!is_element_active_[phase.z[j]])
			{
				active = false;
				break;
			}
		}
		if (active && ((i != P_H2O_LIQUID && i != P_WATER_VAPOR) ||
		               (i == P_H2O_LIQUID && !no_liquid_water) ||
		               (i == P_WATER_VAPOR && no_liquid_water)))
		{
			phase_.push_back(phase);
			Gf_.push_back(phase.model->Gf(phase, T_, P_));
		}
	}
	phase_.shrink_to_fit();
	Gf_.shrink_to_fit();
}

//-----------------------------------------------------------------------------//
   // Construct a state congruent with a Melt_Model and with the given starting phases
State::State(std::string const &name,
             class Melt_Model const &melt,
             unsigned NP,
             unsigned const cphase[],
             double const xphase[])

: name_(name), T_(melt.T()), P_(melt.P()), phase_(melt.phase()), Gf_(melt.Gf())
{
	using namespace std;

	fill(X_, X_+E_END, 0.0);

	ph_[E_H] = P_H2;
	ph_[E_C] = P_GRAPHITE;
	ph_[E_NA] = P_Na;
	ph_[E_MG] = P_Mg;
	ph_[E_AL] = P_Al;
	ph_[E_SI] = P_Si;
	ph_[E_P] = P_P4;
	ph_[E_S] = P_S;
	ph_[E_CL] = P_Cl2;
	ph_[E_K] = P_K;
	ph_[E_CA] = P_Ca;
	ph_[E_TI] = P_Ti;
	ph_[E_CR] = P_Cr;
	ph_[E_MN] = P_Mn;
	ph_[E_FE] = P_Fe;
	ph_[E_ZR] = P_Zr;
	ph_[E_O] = P_O2;

	bool used_[P_END];
	fill(used_, used_+P_END, false);
	for (unsigned i=0; i<NP; ++i)
	{
		Phase const &phase = phase_[cphase[i]];
		unsigned const N = phase.nz;
		for (unsigned j=0; j<N; ++j)
		{
			unsigned z = phase.z[j];
			if (!used_[z])
			{
				ph_[z] = cphase[i];
				X_[z] = xphase[i];
				used_[z] = true;
				break;
			}
		}
	}
}
