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

#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
//#include <fstream>

//#include "kgb/tostring.h"
//#include "kgb/Assert.h"

//#include "config.h"

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"

#include "phase.hh"

//-----------------------------------------------------------------------------//
using namespace std;

unsigned const N = E_END;
	
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
			if (fabs(aGfi)<1.0e-9)
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
			for (unsigned i=0; i<N; ++i)
			{
				Np[i] -= r*left[i];
				if (fabs(Np[i])<1.0e-10)
				{
					if (fabs(left[i])>1.0e-10)
					{
						n = i;
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


template<class Function>
void solve(double const y, double &x, Function f);

static double Gw; 
static double Ga;
static double Rt;
// Paper says in kJ; J much more likely
double const Wnorm = 1.0e-3;
double const Ww_q = 30967.*Wnorm;
double const Ww_cor = -16098.*Wnorm;
double const Ww_fo = 28874.*Wnorm;
double const Ww_fa = 35634.*Wnorm;
double const Ww_wo = 20375.*Wnorm;
double const Ww_sm = -96938.*Wnorm;
double const Ww_kal = 10374.*Wnorm;
double const Wq_cor = -39120.*Wnorm;
double const Wq_fo = 23661.*Wnorm;
double const Wq_fa = 3421.*Wnorm;
double const Wq_wo = -864.*Wnorm;
double const Wq_sm =  -99039.*Wnorm;
double const Wq_kal = -33922.*Wnorm;
double const Wcor_fo = -30509.*Wnorm;
double const Wcor_fa = -32880.*Wnorm;
double const Wcor_wo = -57918.*Wnorm;
double const Wcor_sm = -130785.*Wnorm;
double const Wcor_kal = -25859.*Wnorm;
double const Wfo_fa = -37257.*Wnorm;
double const Wfo_wo = -31732.*Wnorm;
double const Wfo_sm = -41877.*Wnorm;
double const Wfo_kal = 22323.*Wnorm;
double const Wfa_wo = -12917.*Wnorm;
double const Wfa_sm = -90534.*Wnorm;
double const Wfa_kal = 23649.*Wnorm;
double const Wwo_sm = -13247.*Wnorm;
double const Wwo_kal = 17111.*Wnorm;
double const Wsm_kal = 6523.*Wnorm;

unsigned const M_END = 8;  // Current set includes H, Si, Al, Mg, Fe, Ca, Na, K
double W[M_END][M_END] = 
{
	{0.0,      Ww_q,  Ww_cor,    Ww_fo,    Ww_fa,    Ww_wo,   Ww_sm,   Ww_kal},
	{Ww_q,      0.0,  Wq_cor,    Wq_fo,    Wq_fa,    Wq_wo,   Wq_sm,   Wq_kal},
	{Ww_cor, Wq_cor,     0.0,  Wcor_fo,  Wcor_fa,  Wcor_wo, Wcor_sm, Wcor_kal},
    {Ww_fo,   Wq_fo, Wcor_fo,      0.0,   Wfo_fa,   Wfo_wo,  Wfo_sm,  Wfo_kal},
    {Ww_fa,   Wq_fa, Wcor_fa,   Wfo_fa,      0.0,   Wfa_wo,  Wfa_sm,  Wfa_kal},
    {Ww_wo,   Wq_wo, Wcor_wo,   Wfo_wo,   Wfa_wo,      0.0,  Wwo_sm,  Wwo_kal},
    {Ww_sm,   Wq_sm, Wcor_sm,   Wfo_sm,   Wfa_sm,   Wwo_sm,     0.0,  Wsm_kal},
    {Ww_kal, Wq_kal, Wcor_kal, Wfo_kal,  Wfa_kal,  Wwo_kal, Wsm_kal,      0.0},
};

double minimize(double a, double b, double f(double x))
{
	double fa = f(a), fb = f(b);
	double const phi = 0.5*(1+sqrt(5.));

	double c = b - (b-a)/phi;
	double fc =  f(c);
	double d = a + (b-a)/phi;
	double fd = f(d);

	while (fabs(b-a)>1.0e-8)
	{
		if (fc<fd)
		{
			b = d;
			fb = fd;
			d = c;
			fd = fc;
			c = b - (b-a)/phi;
			fc =  f(c);
		}
		else
		{
			a = c;
			fa = fc;
			c = d;
			fc = fd;
			d = a + (b-a)/phi;
			fd = f(d);
		}
	};
	return (fb<fa? b : a);
	
#if 0
	double f0 = f(x0), f1 = f(x1);
	double x2 = 0.5*(x0+x1);
	double f2 = f(x2);

	do 
	{
		// See if the curve is concave up or down
		double det = f0/(x0-x1)/(x0-x2) + f1/(x1-x0)/(x1-x2) + f2/(x2-x0)/(x2-x1);

		if (det > 0.0)
		{
			// concave up; use minimum if possible'
			
		double xm = 0.5*(f0*(x1+x2)/(x0-x1)/(x0-x2) + f1*(x0+x2)/(x1-x0)/(x1-x2) 
		                 + f2*(x0+x1)/(x2-x0)/(x2-x1))/
			(f0/(x0-x1)/(x0-x2) + f1/(x1-x0)/(x1-x2) + f2/(x2-x0)/(x2-x1));

		if (xm <= x0)
		{
			x1 = x2;
			f1 = f2;
			x2 = 0.5*(x0+x1);
			f2 = f(x2);
		}
		else if (xm<=x2)
		{
			x1 = x2;
			f1 = f2;
	      x2 = xm; 
			f2 = f(xm);
		}
		else if (xm<x1)
		{
			x0 = x2;
			f0 = f2;
			x2 = xm;
			f2 = f(xm);
		}
		else // xm >= x1
		{
			x0 = x2;
			f0 = f2;
			x2 = 0.5*(x0+x1);
			f2 = f(x2);
		}
		}
		else
		{
			// concave down; pick lower of endpoints 
			if (f0<f1)
			{
				x1 = x2;
				f1 = f2;
			}
			else
			{
				x0 = x2;
				f0 = f2;
			}
			x2 = 0.5*(x0+x1);
			f2 = f(x2);
		}
	} while (fabs(x1-x0)>1.0e-8);
	return x2;
#endif
}

unsigned NM;
unsigned melt_element[E_END];
double X[M_END];
double p[M_END];
double mGfi[M_END];
double T;

//-----------------------------------------------------------------------------//
double Gfmelt(double const e)
{
	unsigned const NMM = NM-1;
	
	unsigned i0 = melt_element[NMM];
	double x[NM];

	double Gf = 0.0;
	double sum = 0.0;

	for (unsigned i=0; i<NMM; ++i)
	{

		x[i] = max(0.0, min(1.0, X[i] + e*p[i]));
		sum += x[i];
	}
	x[NMM] = max(0.0, min(1.0, 1-sum));
	for (unsigned i=0; i<NM; ++i)
	{
		unsigned const ii = melt_element[i];
		if (x[i]>0.0)
		{
			Gf += mGfi[ii]*x[i] + R*T*x[i]*log(x[i]);
			for (unsigned j=0; j<NM; ++j)
			{
				unsigned const jj = melt_element[j];
				Gf += 0.5*x[i]*x[j]*W[ii][jj];
			}
		}
	}
	return Gf;
}

//-----------------------------------------------------------------------------//
double find_eutectic(double const T)
{
	::T = T;
	  // Compute the initial gradient and Hessian matrix
	unsigned const NMM = NM-1;
	double H[NMM][NMM];
	double grad[NMM];
	unsigned i0 = melt_element[NMM];

	unsigned const N = NMM;
	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_matrix *V = gsl_matrix_alloc(N, N);
	gsl_vector *S = gsl_vector_alloc(N);
	gsl_vector *work = gsl_vector_alloc(N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);	
	gsl_vector *aGf = gsl_vector_alloc(P_END);
	gsl_permutation *permute = gsl_permutation_alloc(P_END);

	DO:

	for (unsigned i=0; i<NMM; ++i)
	{
		unsigned const ii = melt_element[i];

		// Gradient
		
		grad[i] = mGfi[ii] - mGfi[i0] + R*T*(log(X[i]) - log(X[NMM])) + W[i0][ii]*X[NMM];
		for (unsigned j=0; j<NMM; j++)
		{
			unsigned const jj = melt_element[j];
			grad[i] += W[ii][jj]*X[j] - W[i0][jj]*X[j];
		}
		
		// Hessian

		H[i][i] = R*T*(1/X[i] + 1/X[NMM]) - 2*W[i0][ii];
		for (unsigned j=i+1; j<NMM; ++j)
		{
			unsigned const jj = melt_element[j];
			H[i][j] = R*T/X[NMM] - W[i0][ii] + W[ii][jj] - W[i0][jj];
		}
		for (unsigned j=0; j<i; j++)
		{
			H[i][j] = H[j][i];
		}
	}

	// Solve for search direction

	for (unsigned i=0; i<N; ++i)
	{
		for (unsigned j=0; j<N; j++)
		{
			gsl_matrix_set(A, i, j, H[i][j]);
		}
		gsl_vector_set(b, i, -grad[i]);
	}

	gsl_linalg_SV_decomp (A, V, S, work);		
	gsl_linalg_SV_solve (A, V, S, b, x);

	double x0 = -numeric_limits<double>::max();
	double x1 = numeric_limits<double>::max();
	double sum = 0.0;
		bool singular = true;
	for (unsigned i=0; i<NMM; ++i)
	{
		p[i] = gsl_vector_get(x, i);
		if (p[0] != 0.0) singular = false;
		double e = -X[i]/p[i];
		if (e>0) x1 = min(x1, e); else x0 = max(x0, e);
		e = (1-X[i])/p[i];
		if (e>0) x1 = min(x1, e); else x0 = max(x0, e);
		sum -= p[i];
	}
    if (!singular)
	{
		double e = -X[NMM]/sum;
		if (e>0) x1 = min(x1, e); else x0 = max(x0, e);
		e = (1-X[NMM])/sum;
		if (e>0) x1 = min(x1, e); else x0 = max(x0, e);
	}
    else
	{
		// search along inverse grad
		for (unsigned i=0; i<NMM; ++i)
		{
			p[i] = -grad[i];
			if (p[0] != 0.0) singular = false;
			double e = -X[i]/p[i];
			if (e>0) x1 = min(x1, e);
			e = (1-X[i])/p[i];
			if (e>0) x1 = min(x1, e);
			sum -= p[i];
		}
		double e = -X[NMM]/sum;
		if (e>0) x1 = min(x1, e);
		e = (1-X[NMM])/sum;
		if (e>0) x1 = min(x1, e);
	}
		
	double x2 = minimize(x0, x1, Gfmelt);

	if (fabs(x2)>1.0e-4)
	{
		for (unsigned i=0; i<NMM; ++i)
		{
			X[i] += x2*p[i];
		}
		X[NMM] += x2*sum;
		goto DO;
	}
	
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(aGf);
	gsl_permutation_free(permute);
		return Gfmelt(0.0);
}

//-----------------------------------------------------------------------------//
void update(double const T, 
            double const P, 
            vector<Phase> &phase,
            bool const oxygen_specified, 
            bool const oxygen_FMQ,
            double &pO2,
	        double Np[E_END],
	        unsigned p[E_END],
            double volume[N])
{
	vector<double> Gf(P_END);
	initialize_Gf(T, P, phase, Gf);

	// Solve for end point phases.
	double element_activity[E_END];
	do_update(T, P, phase, Gf, oxygen_specified, oxygen_FMQ, pO2, element_activity, Np, p, volume);

	// Now construct activities of melt end members

	unsigned endmember[M_END] = {
		P_WATER_VAPOR,
		P_SiO2_MELT,
		P_CORUNDUM_LIQUID, 
		P_FORSTERITE_LIQUID,
		P_FAYALITE_LIQUID,
		P_WOLLASTONITE_LIQUID,
		P_K_FELDSPAR_LIQUID,
		P_ALBITE_LIQUID};
	
	unsigned endmember_element[M_END] = {
		E_H,
		E_SI,
		E_AL, 
		E_MG,
		E_FE,
		E_CA,
		E_K,
		E_NA};

	for (unsigned i=0; i<M_END; ++i)
	{
		unsigned ie = endmember[i];
		double aGfi = Gf[ie];
		Phase const &phi = phase[ie];
		for (unsigned j=0; j<phi.nz; j++)
		{
			aGfi += phi.n[j]*element_activity[phi.z[j]];
		}
		mGfi[i] = aGfi;
	}

	// Now rebase last two to sodium metasilicate, Na2SiO3, and kalsilite, KAlSi04, as liquids.

	double Asm = 2*mGfi[7] - mGfi[2] - 5*mGfi[1] 
		- (R*T*(2*log(1./7.)/7. + 5*log(5./7.)/7.) + Wcor_sm/49. + Wq_sm*5./49. + Wq_cor*5./49.);
	
	double Akal = 3*mGfi[6] - 2*mGfi[1] 
		- (R*T*(log(1./3.)/3 + 2*log(2./3.)/3) + Wq_kal*2./9.);
	
	mGfi[6] = Akal;
	mGfi[7] = Asm;

	// Now compress out elements not actually present.
    NM = 0;
	for (unsigned i=0; i<E_END; ++i)
	{
		if (Np[i] > 1.e-9)
		{
			// i is an element present in the mix that could go into the melt. Look for it in element melt table.
			for (unsigned j=0; j<M_END; ++j)
			{
				if (endmember_element[j] == i)
				{
					melt_element[NM++] = j;
				}
			}
		}
	}


	// If there's only one liquid, we're done.
	if (NM<2)
		return;
	
	// Minimize.

	// Initial guess is equal quantities of all allowed endmembers.
	fill(X, X+M_END, 0.0);
	for (unsigned i=0; i<NM; ++i)
	{
	  X[i] = 1.0/NM;
	}
	double Geu = find_eutectic(T);

	if (Geu<1.0e-9)
	{
		double x[M_END];
	    fill(x, x+M_END, 0.0);
	    for (unsigned i=0; i<NM; ++i)
	    {
	      x[melt_element[i]] = X[i];
	    }
		double nO = x[0] + 2*x[1] + 3*x[2] + 4*x[3] + 4*x[4] + 3*x[5] + 4*x[6] + 3*x[7];
		double nH = 2*x[0];
		double nSi = x[1] + x[3] + x[4] + x[5] + x[6] + x[7];
		double nAl = 2*x[2] + x[6];
		double nMg = 2*x[3];
		double nFe = 2*x[4];
		double nCa = x[5];
		double nK = x[6];
		double nNa = 2*x[7];
		double W = 16*nO + nH + 28*nSi + 26*nAl + 24*nMg + 56*nFe + 40*nCa + 39*nK + 23*nNa;
		
		double V = phase[P_WATER_VAPOR].V*x[0] 
			+ phase[P_SiO2_MELT].V*x[1]
			+ phase[P_CORUNDUM_LIQUID].V*x[2]
			+ phase[P_FORSTERITE_LIQUID].V*x[3]
			+ phase[P_FAYALITE_LIQUID].V*x[4]
			+ phase[P_WOLLASTONITE_LIQUID].V*x[5]
			+ phase[P_KALSILITE].V*x[6]
		+ 4.68*x[7];

		Geu -= nO*element_activity[E_O];
		Geu -= nH*element_activity[E_H];
		Geu -= nSi*element_activity[E_SI];
		Geu -= nAl*element_activity[E_AL];
		Geu -= nMg*element_activity[E_MG];
		Geu -= nFe*element_activity[E_FE];
		Geu -= nCa*element_activity[E_CA];
		Geu -= nK*element_activity[E_K];
		Geu -= nNa*element_activity[E_NA];

		phase.push_back(Phase{"melt", 9, {E_O, E_H, E_SI, E_AL, E_MG, E_FE, E_CA, E_K, E_NA}, 
			             {nO, nH, nSi, nAl, nMg, nFe, nCa, nK, nNa},
		       	W, V, Geu, 0.0,  0.,0.,0.,0., 0.,0.,0.,0., MAGMA});

		Gf.push_back(Geu);

		do_update(T, P, phase, Gf, oxygen_specified, oxygen_FMQ, pO2, element_activity, Np, p, volume);
	}
	// else we've done all the melting we can
}

