/*
 * melt.cc
 * Copyright (C) 2019 Kent G. Budge <kgb@kgbudge.com>
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
#include <fstream>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"

#include "ds++/Assert.hh"

#include "algorithm.hh"
#include "constants.hh"
#include "element.hh"
#include "melt.hh"
#include "model.hh"
#include "phase.hh"

//-----------------------------------------------------------------------------//
using namespace std;

//-----------------------------------------------------------------------------//
// Paper says in kJ; J much more likely
double const Wnorm = 1.0e-3;
double const Ww_q = 30967.*Wnorm;
double const Ww_cor = -16098.*Wnorm;
double const Ww_fo = 28874.*Wnorm;
double const Ww_fa = 35634.*Wnorm;
double const Ww_wo = 20375.*Wnorm;
double const Ww_sm = -96938.*Wnorm;
double const Ww_kal = 10374.*Wnorm;
double const Wq_cor = -108758*Wnorm;// recalibrated from Q-Cor eutectic at 1346K
double const Wq_fo = 7346*Wnorm; // calibrated from Q-Fo eutectic at 1815K
double const Wq_fa = 11550.*Wnorm; // calibrated from eutectic near 1422K
double const Wq_wo = 14834.*Wnorm;// calibrated from eutectic of 1699K
double const Wq_sm = -6000*Wnorm;// calibrated from eutectic at 1062K
double const Wq_kal = -33922.*Wnorm;
double const Wcor_fo = -30509.*Wnorm;
double const Wcor_fa = -32880.*Wnorm;
double const Wcor_wo = -57918.*Wnorm;
double const Wcor_sm = 100000.*Wnorm; //-130785.*Wnorm;
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

double const W[M_END][M_END] =
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

//-----------------------------------------------------------------------------//
double Gfmelt(double const T,
              unsigned const NM,
              unsigned const melt_element[],
              double const mGfi[],
              double const X[],
              double const p[],
              double const e)
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
			if (ii==M_H2O)
			{
				Gf += R*T*(x[i]*log(x[i]) + (1-x[i])*log(1-x[i]));
			}
		}
	}
	return Gf;
}

//-----------------------------------------------------------------------------//
double find_eutectic(double const T, 
                     unsigned const NM, 
                     unsigned const melt_element[],
                     double const mGfi[],
                     double X[])
{
	  // Compute the initial gradient and Hessian matrix
	unsigned const NMM = NM-1;
	double H[NM][NM];
	double grad[NM];
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

	for (unsigned i=0; i<NM; ++i)
	{
		unsigned const ii = melt_element[i];

		// Gradient

		grad[i] = mGfi[ii] + R*T*(1 + log(X[i]));
		for (unsigned j=0; j<NM; ++j)
		{
			unsigned const jj = melt_element[j];
			grad[i] += 0.5*X[j]*W[ii][jj]; // no special code for i==j since W[ii][jj] is always zero.
		}
		if (ii==M_H2O)
		{
			grad[i] += R*T*(log(X[i]) - log(1-X[i]));
		}
		
		// Hessian

		H[i][i] = R*T/X[i];
		for (unsigned j=0; j<NM; ++j)
		{
			if (i!=j)
			{
				unsigned const jj = melt_element[j];
				H[i][j] = 0.5*W[ii][jj]; 
			}
		}
		if (ii==M_H2O)
		{
			H[i][i] += R*T*(1/X[i] + 1/(1-X[i]));
		}
	}
		
	// Solve for search direction,  imposing constaint X[NMM] = 1 - sum(X[i!=NMM])

	for (unsigned i=0; i<N; ++i)
	{
		for (unsigned j=0; j<N; j++)
		{
			gsl_matrix_set(A, i, j, H[i][j] - 2*H[i][NMM] + H[NMM][NMM]);
		}
		gsl_vector_set(b, i, grad[NMM]-grad[i]);
	}

	gsl_linalg_SV_decomp (A, V, S, work);		
	gsl_linalg_SV_solve (A, V, S, b, x);

	double x0 = -numeric_limits<double>::max();
	double x1 = numeric_limits<double>::max();
	double sum = 0.0;
	bool singular = true;
    double p[M_END];
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

    double x2 = minimize(x0, x1, [=](double const e)
                         {return Gfmelt(T,
                                        NM,
                                        melt_element,
                                        mGfi,
                                        X,
                                        p,
                                        e);});
		
	if (fabs(x2*sum)>1.0e-7)
	{
		sum = 0.0;
		for (unsigned i=0; i<NMM; ++i)
		{
			X[i] += x2*p[i];
			X[i] = min(1.0, max(1e-50, X[i]));
			sum += X[i];
		}
		X[NMM] = min(1.0, max(1e-50, 1-sum));
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

	return Gfmelt(T,
	              NM,
	              melt_element,
	              mGfi,
	              X,
	              p,
	              0.0);
}

//-----------------------------------------------------------------------------//
bool melt(double T, 
            double P, 
            vector<Phase> const &phase,
            vector<double> const &Gf,
            double const element_activity[E_END], 
            double const Np[E_END], 
            unsigned const p[E_END], 
            double &Geu,
            struct Phase &new_phase)
{
	double mGfi[M_END];

	for (unsigned i=0; i<M_END; ++i)
	{
		unsigned ie = endmember[i];
		double aGfi = Gf[ie];
		Phase const &phi = phase[ie];
		for (unsigned j=0; j<phi.nz; j++)
		{
			aGfi += phi.n[j]*element_activity[phi.z[j]];
		}
		mGfi[i] = aGfi*endmember_mole[i];
	}

	// Now compress out elements not actually present.
	bool is_present[E_END];
	fill(is_present, is_present+E_END, false);
	for (unsigned i=0; i<E_END; ++i)
	{
		if (Np[i] > 1.e-9)
		{
			// i is a phase that is present in quantity.
			unsigned j = p[i];
			Phase const &ph = phase[j];
			for (unsigned k=0; k<ph.nz; k++)
			{
				if (ph.n[k]>1.e-9)
				  is_present[ph.z[k]] =  true;
			}
		}
	}
	unsigned NM = 0;
	unsigned melt_element[M_END];
	for (int i=0; i<M_END; ++i)
	{
		if (is_present[endmember_element[i]])
		{
			melt_element[NM++] = i;
		}
	}

	// If there's only one liquid, we're done.
	if (NM<2)
		return false;

	// Minimize.

	// Initial guess is equal quantities of all allowed endmembers.
	double X[M_END];
	fill(X, X+M_END, 0.0);
	for (unsigned i=0; i<NM; ++i)
	{
		X[i] = 1.0/NM;
	}
	Geu = find_eutectic(T, 
	                    NM, 
	                    melt_element,
	                    mGfi,
	                    X);
	
	if (Geu<-1.0e-9)
	{
		double x[E_END];
		fill(x, x+E_END, 0.0);
		double V = 0.0;
		for (unsigned i=0; i<NM; ++i)
		{
			unsigned j = melt_element[i];
			unsigned p = endmember[j];
			Phase const ph = phase[p];
			V += ph.V*X[i]*endmember_mole[j];
			for (unsigned n = 0; n<ph.nz; ++n)
			{
				x[ph.z[n]] += X[i]*ph.n[n]*endmember_mole[j];
			}
		}
		for (unsigned e=0; e<E_END; ++e)
		{
			Geu -= x[e]*element_activity[e];
		}
		new_phase.name = "melt";
		new_phase.nz = M_END+1;
		for (unsigned i=0; i<M_END; ++i)
		{
	      new_phase.z[i] = endmember_element[i];
	      new_phase.n[i] = x[endmember_element[i]];
		}
		new_phase.z[M_END] = E_O;
		new_phase.n[M_END] = x[E_O];
		new_phase.V = V;
		new_phase.model = MELT;
		return true;
	}
	return false;
}

