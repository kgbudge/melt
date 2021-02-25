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

#include <iomanip>
#include <iostream>

#include "Assert.hh"

#include "Model.hh"
#include "element.hh"

//-----------------------------------------------------------------------------//
void 
State::initialize_globals(double T, double P,   
                          double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
                          double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
                          double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
                          double nFe2O3,  double nZrO2)
{
	Require(T>0.0);
	Require(P>0.0);
	
	using namespace std;

	T_ = T;
	P_ = P;

	Require(nH2O>=0.0);
	is_element_active_[E_H] = nH2O>0.0;
	Require(nCO2>=0.0);
	is_element_active_[E_C] = nCO2>0.0;

	is_element_active_[E_O] = nH2O>0.0 || nCO2>0.0 || nNa2O>0.0 || nMgO > 0.0 || nAl2O3 > 0.0 ||
		nSiO2 > 0.0 || nP2O5 > 0.0 || nK2O > 0.0 || nCaO > 0.0 || nTiO2 > 0.0 || nCr2O3 > 0.0 ||
		nMnO > 0.0 || nFeO > 0.0 || nFe2O3 > 0.0 || nZrO2 > 0.0;
	
	Require(nNa2O>=0.0);
	is_element_active_[E_NA] = nNa2O>0.0;
	Require(nMgO>=0.0);
	is_element_active_[E_MG] = nMgO>0.0;
	Require(nAl2O3>=0.0);
	is_element_active_[E_AL] = nAl2O3>0.0;
	Require(nSiO2>=0.0);
	is_element_active_[E_SI] = nSiO2>0.0;
	Require(nP2O5>=0.0);
	is_element_active_[E_P] = nP2O5>0.0;
	Require(nS>=0.0);
	is_element_active_[E_S] = nS>0.0;
	Require(nCl>=0.0);
	is_element_active_[E_CL] = nCl>0.0;
	Require(nK2O>=0.0);
	is_element_active_[E_K] = nK2O>0.0;
	Require(nCaO>=0.0);
	is_element_active_[E_CA] = nCaO>0.0;
	Require(nTiO2>=0.0);
	is_element_active_[E_TI] = nTiO2>0.0;
	Require(nCr2O3>=0.0);
	is_element_active_[E_CR] = nCr2O3>0.0;
	Require(nMnO>=0.0);
	is_element_active_[E_MN] = nMnO>0.0;
	Require(nFeO>=0.0 && nFe2O3>=0.0);
	is_element_active_[E_FE] = nFeO>0.0 || nFe2O3 > 0.0;
	Require(nZrO2>=0.0);
	is_element_active_[E_ZR] = nZrO2>0.0;
	
	NP_ = 0;

	// Basic phases are always active.
	for (unsigned i=0; i<E_END; ++i)
	{
		Phase const &phase = ::phase[i];
		if (i!= phase.index)
		{
			cout << "ERROR: phase index out of synch for " << phase.name << endl;
			exit(1);
		}
		phase_[NP_] = phase;
		Gf_[NP_] = phase.model->Gf(phase, T_, P_);
		NP_++;
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
			phase_[NP_] = phase;
			Gf_[NP_] = phase.model->Gf(phase, T_, P_);
			NP_++;
		}
	}

	// Find all potentiall fusible phases
	 for (unsigned p=0; p<NP_; ++p)
	 {
		 Phase const &ph = phase_[p];	
		 cout << ph.name << endl;
		 unsigned const N = ph.nz; // Number of elements in the phase
		 double xO = 0.0;          // To accumulate oxygen balance of the phase
		 double xFe = 0.0;
		 bool fusible = true;
		 for (unsigned j=0; j<N && fusible; ++j) 
		 {
			 double const xj = ph.n[j];  // Number of moles of element j in the phase.
			 unsigned const z = ph.z[j];
			 fusible = fusible && (is_element_active_[z]);
			 switch(z)             // Switch on the element atomic number
			 {
				 case E_H:
					 xO -= 0.5*xj;
					 break; 

				 case E_C:
					 xO -= 2*xj;
					 break;

				 case E_O:
					 xO += xj;
					 break;

				 case E_NA:
					 xO -= 0.5*xj;
					 break;

				 case E_MG:
					 xO -= xj;
					 break;

				 case E_AL:
					 xO -= 1.5*xj;
					 break;

				 case E_SI:
					 xO -= 2*xj;
					 break;

				 case E_S:
					 break;

				 case E_CL:
					 xO += 0.5*xj;
					 break;

				 case E_K:
					 xO -= 0.5*xj;
					 break;

				 case E_CA:
					 xO -= xj;
					 break;

				 case E_FE:
					 xFe += xj;
					 break;

				 default:
					 // non-fusible mineral
					 cout << "phase " << ph.name << " cannot melt." << endl;
					 fusible = false;
					 break;
			 }
		 }
	     is_fusible_[p] = fusible && fabs(xO-xFe)<1e-9;
		 // at present, cannot handle ferric or oxidized sulfur melts
	 }

	phase_[NP_] = phase_[0];
	Gf_[NP_] = 1.0e10;
	is_fusible_[NP_] = false;
	NP_++;

	 // Compute end-member melt phase free energies
	 for (unsigned i=0; i<M_END; ++i)
	 {
		 Phase const &ph = ::phase[melt_endmember[i]];
		 Gfm_[i] = ph.model->Gf(ph, T_, P_);
	 }
	
	gsl_matrix_free(gsl_A_);
	gsl_matrix_free(gsl_V_);
	gsl_vector_free(gsl_S_);
	gsl_vector_free(gsl_work_);
	gsl_vector_free(gsl_b_);
	gsl_vector_free(gsl_x_);
	gsl_vector_free(gsl_aGf_);
	gsl_permutation_free(gsl_permute_);
	
	gsl_A_ = gsl_matrix_alloc(E_END, E_END);
	gsl_V_ = gsl_matrix_alloc(E_END, E_END);
	gsl_S_ = gsl_vector_alloc(E_END);
	gsl_work_ = gsl_vector_alloc(E_END);
	gsl_b_ = gsl_vector_alloc(E_END);
	gsl_x_ = gsl_vector_alloc(E_END);	
	gsl_aGf_ = gsl_vector_alloc(NP_);
	gsl_permute_ = gsl_permutation_alloc(NP_);
}
