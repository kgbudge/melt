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
#include "model.hh"

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
	cout << "Potentially active phases:" << endl;
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
			cout << phase[i].name << " " << fixed << setprecision(3) << (1000*Gf[i]/amu[i]) << " kj/kg" << endl;
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
void do_update(double const T, 
               double const P, 
               vector<Phase> const &phase,
               vector<double> Gf,
               bool const oxygen_specified, 
               bool const oxygen_FMQ,
               double &pO2,
               double element_activity[E_END],
               double Np[E_END],
               unsigned p[E_END],
               double volume[N])
{
	unsigned const P_END = phase.size();
	double GfO2;
	if (oxygen_specified)
	{
		if (oxygen_FMQ)
		{
			GfO2 = 2*Gf[P_MAGNETITE]+3*Gf[P_SiO2_QUARTZ]-3*Gf[P_FAYALITE];
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
	
	// now iterate to find composition

	for(;;)
	{
		// find absolute free energy of elements

		// Build a linear system and solve for element free energies
		
		gsl_matrix_set_zero(A);
		gsl_vector_set_zero(b);

		if (oxygen_specified)
		{
			Gf[P_O2] = GfO2;
		}

		for (unsigned i=0; i<N; ++i)
		{
			unsigned const pi = p[i];
			gsl_vector_set(b, i, -Gf[pi]);
			unsigned const nz = phase[pi].nz;
			for (unsigned j=0; j<nz; j++)
			{
				gsl_matrix_set(A, i, phase[pi].z[j], phase[pi].n[j]);
			}
		}

		gsl_linalg_SV_decomp (A, V, S, work);		
		gsl_linalg_SV_solve (A, V, S, b, x);

		// Compute absolute free energies of phases

		for (unsigned i=0; i<E_END; ++i)
		{
			element_activity[i] = gsl_vector_get(x, i);
		}

		for (unsigned i=0; i<P_END; ++i)
		{
			double aGfi = Gf[i];
			for (unsigned j=0; j<phase[i].nz; j++)
			{
				aGfi += phase[i].n[j]*element_activity[phase[i].z[j]];
			}
			gsl_vector_set(aGf, i, aGfi);
		}
	
		gsl_sort_vector_index(permute, aGf);

		// Find a composition basis for the current composition.

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
		
		for (unsigned ii=0; ii<P_END; ++ii)
		{
			unsigned const ip = gsl_permutation_get(permute, ii);
			if (oxygen_specified && ip == P_O2)
				continue;

			double aGfi = gsl_vector_get(aGf, ip); 
			if (fabs(aGfi)<1.0e-9 || !isfinite(aGfi))
			{
				goto DONE;
			}

			// Construct an equation that reduces the Gibbs function further.

			// Construct the reaction

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

			// Can the reaction proceed?

			double r = numeric_limits<double>::max();
			double rGf = Gf[ip];
			for (unsigned i=0; i<N; ++i)
			{
				if (left[i]>1.0e-10)
				{
					r = min(r, Np[i]/left[i]);
				}
				rGf -= left[i]*Gf[p[i]];
			}

			if (r<=1.0e-10)
			{
				continue;
			}

			// Yes, the reaction has somewhere to go

			cout << "Performing reaction ";
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
					cout << phase[p[i]].name;
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
					cout << phase[p[i]].name;
				}
			}
			cout << endl;

			unsigned n;
			double npm = 0.0;
			for (unsigned i=0; i<N; ++i)
			{
				Np[i] -= r*left[i];
				if (fabs(Np[i])<1.0e-10)
				{
					if (r*left[i]>npm)
					{
						n = i;
						npm = r*left[i];
					}
					Np[i] = 0.0;
				}
			}
			Np[n] = r;	
			p[n] = ip;

			break;
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
		pO2 = phase[P_O2].model->P(phase[P_O2], T, -2*element_activity[E_O]);
	}
		
	// Convert to volume fraction

	double Vtot = 0.0;
	for (unsigned i=0; i<N; ++i)
	{
		unsigned const pi = p[i];
        volume[i] = Np[i]*phase[pi].model->volume(phase[pi], T, P);
		Vtot += volume[i];
	}
		
	double rnorm = 100.0/Vtot;
	for (unsigned i=0; i<N; ++i)
	{
		volume[i] *= rnorm;
	}
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

	bool is_element_active[E_END];
	fill(is_element_active, is_element_active+E_END, false);
	for (unsigned e=0; e<E_END; ++e)
	{
		if (state.x[e]>0.0)
		{
		    Phase const &ph = phase[state.p[e]];
			for (unsigned i=0; i<ph.nz; ++i)
			{
				is_element_active[ph.z[i]] = true;
			}
		}
	}
	cout << "Active elemental phases:" << endl;
	for (unsigned e=0; e<E_END; ++e)
	{
		if (is_element_active[e])
		{
			cout << " " << phase[state.p[e]].name << endl;
		}
	}
	
	vector<double> Gf, amu;
	calculate_Gf(T, P, is_element_active, phase, Gf, amu);

	double Gftot = 0.0;
	double Mtot = 0.0;
	for (unsigned e=0; e<E_END; ++e)
	{
		Gftot += state.x[e]*Gf[state.p[e]];
		Mtot += state.x[e]*amu[state.p[e]];
	}
	cout << "Initial free energy of formation = " << (1000*Gftot/Mtot) << " kJ/kg" << endl;

	// Solve for end point phases.
	double element_activity[E_END];

	cout << endl << "Ladder search" << endl;
	do_update(T, 
	          P, 
	          phase, 
	          Gf, 
	          oxygen_specified, 
	          oxygen_FMQ, 
	          pO2, 
	          element_activity, 
	          state.x, 
	          state.p, 
	          state.V);
	
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
	cout << "Free energy of formation = " << fixed << setprecision(3) << (1000*Gftot/Mtot) << " kJ/kg" << endl;

	for (unsigned i=0; i<30; ++i)
	{
		Phase new_phase;
		double Geu;
		if (melt(T, P, phase, Gf, element_activity, state.x, state.p, Geu, new_phase))
		{
			Gf.push_back(Geu);
			phase.push_back(new_phase);

			do_update(T, 
			          P, 
			          phase, 
			          Gf, 
			          oxygen_specified, 
			          oxygen_FMQ, 
			          pO2, 
			          element_activity,
			          state.x, 
			          state.p, 
			          state.V);
		}
		// else we've done all the melting we can
		else
		{
			return;
		}
	}
}

