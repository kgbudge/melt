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
             double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
             double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
             double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
             double nFe2O3,  double nZrO2)

: name_(name)
{
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
}

//-----------------------------------------------------------------------------//
State::State(std::string const &name, 
             double const x0[E_END])

: name_(name)
{
	Require(std::count_if(x0, x0+E_END, [](double x){return x<0.0;}) == 0);
	
	using namespace std;

	copy(x0, x0+E_END, X_);

	ph_[E_H] = P_H2;
	X_[E_H] *= 0.5;  // to H2
	ph_[E_C] = P_GRAPHITE;
	ph_[E_NA] = P_Na;
	ph_[E_MG] = P_Mg;
	ph_[E_AL] = P_Al;
	ph_[E_SI] = P_Si;
	ph_[E_P] = P_P4;
	X_[E_P] *= 0.25;  // to P4
	ph_[E_S] = P_S;
	ph_[E_CL] = P_Cl2;
	X_[E_CL] *= 0.5;   // to Cl2
	ph_[E_K] = P_K;
	ph_[E_CA] = P_Ca;
	ph_[E_TI] = P_Ti;
	ph_[E_CR] = P_Cr;
	ph_[E_MN] = P_Mn;
	ph_[E_FE] = P_Fe;
	ph_[E_ZR] = P_Zr;
	ph_[E_O] = P_O2;
	X_[E_O] *= 0.5;   // to O2
}
