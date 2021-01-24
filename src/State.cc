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

#include "phase_enum.hh"

//-----------------------------------------------------------------------------//
State::State(std::string const &name, 
             double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
             double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
             double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
             double nFe2O3,  double nZrO2)

: name(name)
{
	using namespace std;

	p[E_H] = P_H2;
	x[E_H] = nH2O;
	double nO2 = 0.5*nH2O;

	p[E_C] = P_GRAPHITE;
	x[E_C] = nCO2;
	nO2 += nCO2;

	p[E_NA] = P_Na;
	x[E_NA] = 2*nNa2O;
	nO2 += 0.5*nNa2O;

	p[E_MG] = P_Mg;
	x[E_MG] = nMgO;
	nO2 += 0.5*nMgO;

	p[E_AL] = P_Al;
	x[E_AL] = 2*nAl2O3;
	nO2 += 1.5*nAl2O3;

	p[E_SI] = P_Si;
	x[E_SI] = nSiO2;
	nO2 += nSiO2;

	p[E_P] = P_P4;
	x[E_P] = 0.5*nP2O5;
	nO2 += 2.5*nP2O5;

	p[E_S] = P_S;
	x[E_S] = nS;

	p[E_CL] = P_Cl2;
	x[E_CL] = 0.5*nCl;

	p[E_K] = P_K;
	x[E_K] = 2*nK2O;
	nO2 += 0.5*nK2O;

	p[E_CA] = P_Ca;
	x[E_CA] = nCaO;
	nO2 += 0.5*nCaO;

	p[E_TI] = P_Ti;
	x[E_TI] = nTiO2;
	nO2 += nTiO2;

	p[E_CR] = P_Cr;
	x[E_CR] = 2*nCr2O3;
	nO2 += 1.5*nCr2O3;

	p[E_MN] = P_Mn;
	x[E_MN] = nMnO;
	nO2 += 0.5*nMnO;

	p[E_FE] = P_Fe;
	x[E_FE] = nFeO + 2*nFe2O3;
	nO2 += 0.5*nFeO + 1.5*nFe2O3;

	p[E_ZR] = P_Zr;
	x[E_ZR] = nZrO2;
	nO2 += nZrO2;

	p[E_O] = P_O2;
	x[E_O] = nO2;
}
