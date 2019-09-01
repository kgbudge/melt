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
#include "D2.hh"
#include "element.hh"
#include "melt.hh"
#include "model.hh"
#include "phase.hh"
#include "State.hh"

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
double const Wq_cor = 33470*Wnorm; //1862K 32423*Wnorm;//1839K 28409*Wnorm; calibrated from SiO2-Al2O3 eutectic at 1868K
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
class Melt_Model
{
	public:
		Melt_Model(double const T, 
		           double const P, 
		           vector<Phase> const &phase,
		           vector<double> const &Gf,
		           State const &state) noexcept(false);

		unsigned nm() const noexcept { return nm_; }
		double x(unsigned i) const { Require(i<nm()); return x_[i]; }

		template<typename Real>
		Real Gfmelt(vector<double> const &x) const;

		template<typename Real>
		Real Gfm(vector<Real> const &x) const;

		double Gfmelt(std::vector<double> const &X, double p[M_END], double e) const;
		
		Phase minimize_Gf(vector<double> &x);

		

	private:
		double T_;
		unsigned nm_;
		double basis_[E_END][M_END];
		unsigned p_[M_END];
		double x_[M_END];
		double Gfs_[E_END];
		double Gfm_[M_END];
		double Gf0_; // non-meltable phases total free energy
};

//-----------------------------------------------------------------------------//
/*! Create a Melt_Model reflecting possible melting of an existing State.
 *
 * \param T Temperature (K)
 * \param P Pressure (kbar)
 * \param phase Current phase library
 * \param Gf Gibbs free energy of each phase at T and P (kJ/mol)
 * \param state Current state of the mineral ensemble.
 */ 
 Melt_Model::Melt_Model(double const T, 
                       double const P, 
                       vector<Phase> const &phase,
                       vector<double> const &Gf,
                       State const &state)
{
	// First calculate degrees of freedom. Each active phase that *can* melt
	// is a degree of freedom, allowed to vary from 0 (no melting) to the quantity
	// of that phase. We calculate the amount of each basic melt phase produced
	// by a mole of each meltable phase. Any non-fusible mineral has its free energy
	// added to Gf0_ and is compressed out of the fusible phase set.

	T_ = T;  
	nm_ = 0;     // Number of fusible phases
	Gf0_ = 0.0;  // Gibbs free energy contribution of nonfusible phases 
	for (unsigned i=0; i<E_END; ++i)
	{
		if (state.x[i]>0.0)  // Is this phase actually present?
		{
			fill(basis_[nm_], basis_[nm_]+M_END, 0.);
			unsigned p = state.p[i];  // Prepare to compute a melt basis for the candidate phase.
			p_[nm_] = p;              // Save a candidate fusible phase phase index
			x_[nm_] = state.x[i];     // Save a candidate fusible phase quantity
			Gfs_[nm_] = Gf[p];        // Save the candidate phase unmelted Gibbs free energy
			Phase const &ph = phase[p];	
			unsigned const N = ph.nz; // Number of elements in the phase
			double xO = 0.0;          // Oxygen balance of the basis of the phase
			bool really_solid = false; // Initial assumption is that the phase is fusible
			for (unsigned j=0; j<N && !really_solid; ++j) // Try to build a basis, element by element.
			{
				double const xj = ph.n[j];  // Number of moles of element j in the phase.
				switch(ph.z[j])             // Switch on the element atomic number to build basis
				{
					case E_H:
						xO -= 0.5*xj;
						basis_[nm_][M_H2O] += 0.5*xj;
						break; 

				    case E_C:
						xO -= 2*xj;
						basis_[nm_][M_CO2] += xj;
						break;
						
					case E_O:
						xO += xj;
						break;

					case E_NA:
						xO -= 0.5*xj;
						basis_[nm_][M_Na2SiO3] += 0.5*xj;
						basis_[nm_][M_SiO2] -= 0.5*xj;
						break;

					case E_MG:
						xO -= xj;
						basis_[nm_][M_MgO] += 0.5*xj;
						break;

					case E_AL:
						basis_[nm_][M_Al2O3] += 0.5*xj;
						xO -= 1.5*xj;
						break;

					case E_SI:
						basis_[nm_][M_SiO2] += xj;
						xO -= 2*xj;
						break;

					case E_S:
						basis_[nm_][M_S2] += xj;
						break;

					case E_CL:
						basis_[nm_][M_NaCl] += xj;
						basis_[nm_][M_Na2SiO3] -= 0.5*xj;
						basis_[nm_][M_SiO2] += 0.5*xj;
						xO += 0.5*xj;
						break;

					case E_K:
						basis_[nm_][M_KAlSi2O6] += xj;
						xO -= 0.5*xj;
						basis_[nm_][M_Al2O3] -= 0.5*xj;
						basis_[nm_][M_SiO2] -= 2*xj;
						break;

					case E_CA:
						xO -= xj;
						basis_[nm_][M_CaO] += xj;
						break;

					case E_FE:
						xO -= xj;
						basis_[nm_][M_Fe2SiO4] += 0.5*xj;
						break;

					default:
						// non-fusible mineral
						cout << "phase " << ph.name << " cannot melt." << endl;
						really_solid = true;
						break;
				}
			}
			if (really_solid || fabs(xO)>1e-9) // at present, cannot handle ferric or oxidized sulfur melts
			{
				Gf0_ += state.x[i]*Gf[p];
			}
			else
			{
				nm_++;  // Accept this candidate phase; it's fusible.
			}
		}
	}
	for (unsigned i=0; i<M_END; ++i)
	{
		Gfm_[i] = Gf[endmember[i]];    // Save the Gibbs free energy of each melt endmember.
	}
}

//-----------------------------------------------------------------------------//
/*! Calculate molar Gibbs free energy of melt.
 * 
 * \param x Mole fraction of each melt member. Energy returned will be for the
 * amount in x.
 */
template<typename Real>
Real Melt_Model::Gfm(std::vector<Real> const &X) const
{
	unsigned const NM = nm_;

	Real zero = to_Real<Real>(0.0, NM);

	std::vector<Real> x = X;

	// Reorganize in analogy to CIPW norm into preferred melt phases.

	// I have no calcite melt

	// Leucite to orthoclase
    Real Q = x[M_KAlSi2O6];
	x[M_KAlSi3O8] = Q;
	x[M_KAlSi2O6] = zero;
	x[M_SiO2] -= Q;

	// Sodium metasilicate to albite
	Q = min(x[M_Na2SiO3], x[M_Al2O3]);
	x[M_NaAlSi3O8] += 2*Q;
	x[M_Na2SiO3] -= Q;
	x[M_Al2O3] -= Q;
	x[M_SiO2] -= 5*Q;

	// Anorthite
	Q = min(x[M_CaO], x[M_Al2O3]);
	x[M_CaAl2Si2O8] = Q;
	x[M_CaO] -= Q;
	x[M_Al2O3] -= Q;
	x[M_SiO2] -= 2*Q;

	// Acmite should come next, but I have no melt for it.

	// Oxidized iron (magnetite and hematite) should come next but I have no melts for them.

	// Magnesium diopside. I have no hedenbergite melt.

	Q = min(x[M_CaO], x[M_MgO]);
	x[M_CaMgSi2O6] = Q;
	x[M_CaO] -= Q;
	x[M_MgO] -= Q;
	x[M_SiO2] -= 2*Q;

	// Wollastonite

	Q = x[M_CaO];
	x[M_CaSiO3] = Q;
	x[M_CaO] = zero;
	x[M_SiO2] -= Q;

	// Magnesium pyroxene (enstatite). I have no ferrosilite melt.

	Q = 2*x[M_MgO];
	x[M_Mg2Si2O6] = Q;
	x[M_MgO] = zero;
	x[M_SiO2] -= 2*Q;

	if (x[M_SiO2]<=0.0)
	{
		// Silica undersaturated

		Real D = -x[M_SiO2];
		x[M_SiO2] = zero;

		// Enstatite to olivine.
		Q = min(D, x[M_Mg2Si2O6]);
		x[M_Mg2SiO4] = Q;
		x[M_Mg2Si2O6] -= Q;
		D -= Q;

		if (D>0.0)
		{
			// Albite to nepheline

			Q = min(0.5*D, x[M_NaAlSi3O8]);
			x[M_NaAlSiO4] += Q;
			x[M_NaAlSi3O8] -= Q;
			D -= 2*Q;

			if (D>0.0)
			{
				// Orthoclase to leucite

				Q = min(x[M_KAlSi3O8], D);
				x[M_KAlSi2O6] += Q;
				x[M_KAlSi3O8] -= Q;
				D -= Q;

				if (D>0.0)
				{
					// Do not have calcium orthosilicate. Go to lime instead.

					// Diopside to wollastonite and olivine.

					Q = min(x[M_CaMgSi2O6], 2*D);
					x[M_CaMgSi2O6] -= Q;
					x[M_CaSiO3] += Q;
					x[M_Mg2SiO4] += 0.5*Q;
					D -= 0.5*Q;

					if (D>0.0)
					{
						// Enstatite to olivine.

						Q = min(D, x[M_Mg2Si2O6]);
						x[M_Mg2Si2O6] -= Q;
						x[M_Mg2SiO4] += Q;
						D -= Q;

						if (D>0.0)
						{
							// Wollastonite to lime

							Q = min(D, x[M_CaSiO3]);
							x[M_CaO] += Q;
							x[M_CaSiO3] -= Q;
							D -= Q;

							if (D>0.0)
							{
								// Olivine to periclase

								Q = min(2*D, 2*x[M_Mg2SiO4]);
								x[M_MgO] += Q;
								x[M_Mg2SiO4] -= 0.5*Q;
								D -= 0.5*Q;
								
								// That's all I have melts for. If quartz is still deficient, this is not a valid melt.
								if (D<-1.0e-9)
								{
									return to_Real<Real>(1e10, NM); // huge Gf turns off this melt.
								}
							}
						}
					}
				}
			}
		}
	}
	
	Real sum = zero;
	for (unsigned i=0; i<M_END; ++i)
	{
		sum += x[i];
	}
	Real rsum = 1.0/(sum + std::numeric_limits<double>::min());
	for (unsigned i=0; i<M_END; ++i)
	{
		x[i] *= rsum;
	}

	// Compute melt free energy.

	Real Gfm = zero;
	for (unsigned i=0; i<M_END; ++i)
	{
		if (x[i]>0.0)
		{
			Gfm += Gfm_[i]*x[i] + R*T_*x[i]*log(x[i]);
			for (unsigned j=0; j<M_END; ++j)
			{
				if (x[j]>0.0)
				{
					Gfm += 0.5*x[i]*x[j]*W[i][j];
				}
			}
/*			if (ii==M_H2O)
			{
				Gfm += R*T_*(x[i]*log(x[i]) + (1-x[i])*log(1-x[i]));
			}
*/		}
	}

	// Now compute total free energy. 
	return sum*Gfm;
}

//-----------------------------------------------------------------------------//
template<typename Real>
Real Melt_Model::Gfmelt(std::vector<double> const &X) const
{
	unsigned const NM = nm_;

	Real zero = to_Real<Real>(0.0, NM);
	std::vector<Real> x(M_END, zero);
	Real Gf = to_Real<Real>(Gf0_, NM);

	for (unsigned i=0; i<NM; ++i)
	{
		Real const xi = to_Real<Real>(min(x_[i], max(0.0, X[i])), i, NM);
		for (unsigned j=0; j<M_END; ++j)
		{
			x[j] += xi*basis_[i][j];
		}
		Gf += (x_[i]-xi)*Gfs_[i];
	}

	// Now compute total free energy. 
	return Gf + Gfm(x);
}

//-----------------------------------------------------------------------------//
/*! Find the degree of melting of each solid phase that minimized total Gibbs free energy.
 * 
 * \param[in,out] X Contains initial guess of amount of melt of each phase. On return,
 * contains amount of melt of each phase that minimizes free energy.
 */
Phase Melt_Model::minimize_Gf(vector<double> &X)
{
	  // Compute the initial gradient and Hessian matrix
	unsigned const N = nm_;

	double H[N][N];
	double grad[N];

	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_matrix *V = gsl_matrix_alloc(N, N);
	gsl_vector *S = gsl_vector_alloc(N);
	gsl_vector *work = gsl_vector_alloc(N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);	
	gsl_vector *aGf = gsl_vector_alloc(P_END);
	gsl_permutation *permute = gsl_permutation_alloc(P_END);

	DO:

		// Calculate initial Gf, gradient, and Hessian.

	D2 Gf2 = Gfmelt<D2>(X);

#if 0 //D2 test
	auto Xp = X;
	Xp[0] -= 0.01;
	Xp[1] -= 0.01;
	D2 Gfmm = Gfmelt<D2>(Xp);
    Xp = X;
	// Xp[0] unchanged
	Xp[1] -= 0.01;
	D2 Gf0m = Gfmelt<D2>(Xp);
	Xp = X;
	Xp[0] += 0.01;
	Xp[1] -= 0.01;
	D2 Gfpm = Gfmelt<D2>(Xp);

	Xp = X;
	Xp[0] -= 0.01;
	// Xp[1] unchanged
	D2 Gfm0 = Gfmelt<D2>(Xp);
    Xp = X;
	// Xp[0] unchanged
	// Xp[1] unchanged
	D2 Gf00 = Gfmelt<D2>(Xp);
	Xp = X;
	Xp[0] += 0.01;
	// Xp[1] unchanged
	D2 Gfp0 = Gfmelt<D2>(Xp);

	Xp = X;
	Xp[0] -= 0.01;
	Xp[1] += 0.01;
	D2 Gfmp = Gfmelt<D2>(Xp);
    Xp = X;
	// Xp[0] unchanged
	Xp[1] += 0.01;
	D2 Gf0p = Gfmelt<D2>(Xp);
	Xp = X;
	Xp[0] += 0.01;
	Xp[1] += 0.01;
	D2 Gfpp = Gfmelt<D2>(Xp);

	double d1 = (Gfp0.y() - Gfm0.y())*50.;
	double d1a = Gf00.dydx(0);
	double d2 = (Gf0p.y() - Gf0m.y())*50.;
	double d2a = Gf00.dydx(1);

	double d11 = (Gfp0.y() - 2*Gf00.y() + Gfm0.y())*10000.;
	double d11a = Gf00.d2ydx2(0, 0);
	double d22 = (Gf0p.y() - 2*Gf00.y() + Gf0m.y())*10000.;
	double d22a = Gf00.d2ydx2(1, 1);
	double d12 = ((Gfpp.y() - Gfmp.y())*50. - (Gfpm.y() - Gfmm.y())*50.)*50.;
	double d12a = Gf00.d2ydx2(0, 1);
#endif

    double Gf = Gf2.y();

	cout << "  Melt parameters:" << endl;
	for (unsigned i=0; i<N; ++i)
	{
		cout << "  X[" << i << "] = " << X[i] << endl;
		
		// Gradient

		grad[i] = Gf2.dydx(i);

		// Hessian

		for (unsigned j=0; j<N; ++j)
		{
			H[i][j] = Gf2.d2ydx2(i, j);
		}
	}
	cout << "  Gf = " << Gf << endl;
		
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

		// Prune singular values

	bool singular = true;
    for (unsigned i=0; i<N; ++i)
	{
		double w = gsl_vector_get(S, i);
		if (fabs(w)<1.0e-9)
		{
			gsl_vector_set(S, i, 0.0);
		}
		else
		{
			singular = false;
		}
	}
		
	gsl_linalg_SV_solve (A, V, S, b, x);

	double x0 = -numeric_limits<double>::max();
	double x1 = numeric_limits<double>::max();
    double p[M_END];
	double norm = 0.0;
	double cnorm = 0.0;
	for (unsigned i=0; i<N; ++i)
	{
		p[i] = gsl_vector_get(x, i);
		// Prune value up against a compositional constraint
		if (p[i]>0.0 && fabs(X[i]-x_[i])<1.0e-9)
		{
			X[i] = x_[i];
			p[i] = 0.0;
		}
		if (p[i]<0.0 && fabs(X[i])<1.0e-9)
		{
			X[i] = 0.0;
			p[i] = 0.0;
		}
		double e = -X[i]/p[i];
		if (p[i]<0) x1 = min(x1, e); else x0 = max(x0, e);
		e = (x_[i]-X[i])/p[i];
		if (p[i]>0) x1 = min(x1, e); else x0 = max(x0, e);
		norm += p[i]*p[i];
		cnorm += x_[i]*x_[i];
	}
	norm = sqrt(norm);
	cnorm = sqrt(cnorm);

	if (!singular && norm > 1.0e-7*cnorm)
	{
		double x2 = minimize(x0, x1, [&](double const e)
			                     {return Gfmelt(X,
			                                    p,
			                                    e);});

		if (fabs(x2)*norm>1.0e-7*cnorm)
		{
			double sum = 0.0;
			for (unsigned i=0; i<N; ++i)
			{
				X[i] += x2*p[i];
				X[i] = min(x_[i], max(0.0, X[i]));
				sum += x2*p[i];
			}
			goto DO;
		}
	}
		
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_vector_free(aGf);
	gsl_permutation_free(permute);

	{
		unsigned const NM = nm_;

		std::vector<double> x(M_END, 0.0);

		// Calculate composition of melt in terms of end members

		for (unsigned i=0; i<NM; ++i)
		{
			double const xi = min(x_[i], max(0.0, X[i]));
			for (unsigned j=0; j<M_END; ++j)
			{
				x[j] += xi*basis_[i][j];
			}
		}

		double sum = 0.0;
		for (unsigned i=0; i<M_END; ++i)
		{
			sum += x[i];
		}
		double rsum = 1.0/(sum + std::numeric_limits<double>::min());
		for (unsigned i=0; i<M_END; ++i)
		{
			x[i] *= rsum;
		}

		// Convert end member composition to elemental composition.

		double xe[E_END];
		fill(xe, xe+E_END, 0.0);
		double V = 0.0;
		for (unsigned i=0; i<M_END; ++i)
		{
			double xi = x[i];
			if (xi>0.0)
			{
				Phase const &ph = phase[endmember[i]];
				unsigned nz = ph.nz;
				for (unsigned j=0; j<nz; ++j)
				{
					xe[ph.z[j]] += xi*ph.n[j];
				}
			    V += ph.V*xi*endmember_mole[i];
			}
		}

		double Gfp = Gfm(x);
		
		Phase Result;
		Result.index = 0;
		Result.name = "mixed melt";
		Result.nz = 0;
		Result.V = V;
		for (unsigned i=0; i<E_END; ++i)
		{
			if (xe[i]>0.0)
			{
				Result.z[Result.nz] = i;
				Result.n[Result.nz] = xe[i];
				Result.nz++;
			}
		}

		Result.S0 = 0.0;
		Result.Hf0 = Gfp;
		Result.model = MELT;

		return Result;
	}
}

//-----------------------------------------------------------------------------//
double Melt_Model::Gfmelt(std::vector<double> const &X, double p[M_END], double e) const
{
	unsigned const N = nm_;
	std::vector<double> x(N);
	for (unsigned i=0; i<N; ++i)
	{
		x[i] = X[i] + e*p[i];
	}
	return Gfmelt<double>(x);
}

//-----------------------------------------------------------------------------//
double melt(double const T, 
            double const P, 
            vector<Phase> const &phase,
            vector<double> const &Gf,
            State const &state, 
            struct Phase &new_phase)
{
	Melt_Model model(T, P, phase, Gf, state);

    // Initial guess is that all phases are half melted. We may try a number of initial guesses eventually.
	unsigned const NM = model.nm();
	vector<double> X(NM);
	for (unsigned i=0; i<NM; ++i)
	{
		X[i] = 0.5*model.x(i);
	}

	new_phase = model.minimize_Gf(X);

	double Gff = model.Gfmelt<double>(X);

	return Gff;
}

