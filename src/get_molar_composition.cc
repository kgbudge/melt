/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * get_molar_composition.cc
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


#include "kgb/tostring.h"

#include "gui.hh"

//-----------------------------------------------------------------------------//
double 
get_molar_composition(std::string &name,
	double &nSiO2,
	double &nTiO2,
	double &nAl2O3,
	double &nFe2O3,
	double &nFeO,
	double &nMnO,
	double &nMgO,
	double &nCaO,
	double &nNa2O,
	double &nK2O,
	double &nP2O5,
	double &nS,
	double &nCr2O3,
	double &nZrO2,
	double &nH2O,
	double &nCO2,
	double &nCl)
{
	using namespace std;

	bool is_molar;

	double total = get_composition(name,
	                               is_molar,
	                               nSiO2,
	                               nTiO2,
	                               nAl2O3,
	                               nFe2O3,
	                               nFeO,
	                               nMnO,
	                               nMgO,
	                               nCaO,
	                               nNa2O,
	                               nK2O,
	                               nP2O5,
	                               nS,
	                               nCr2O3,
	                               nZrO2,
	                               nH2O,
	                               nCO2,
	                               nCl);

	if (total==0.0) return 0.0;
		
	// Calculate molar composition if composition is by weight

	if (button_byweight->get_active())
	{
		nSiO2 /= 60.085;
		nTiO2 /= 79.899;
		nAl2O3 /= 101.961;
		nFe2O3 /= 159.692;
		nFeO /= 71.85;
		nMnO /= 70.94;
		nMgO /= 40.311;
		nCaO /= 56.08;
		nNa2O /= 61.98;
		nK2O /= 94.20;
		nP2O5 /= 109.95;
		nS /= 32.06;
		nCr2O3 /= 223.84;
		nZrO2 /= 123.22;
		nH2O /= 18.02;
		nCO2 /= 44.01;
		nCl /= 35.45;
	}
	
	total = nSiO2 + nTiO2 + nAl2O3 + nFe2O3 + nFeO + nMnO + nMgO + nCaO +
		nNa2O + nK2O + nP2O5 + nS + nCr2O3 + nZrO2 + nH2O + nCO2 + nCl;

	double rnorm = 100.0/total;
	nSiO2 *= rnorm;
	nTiO2 *= rnorm;
	nAl2O3 *= rnorm;
	nFe2O3 *= rnorm;
	nFeO *= rnorm;
	nMnO *= rnorm;
	nMgO *= rnorm;
	nCaO *= rnorm;
	nNa2O *= rnorm;
	nK2O *= rnorm;
	nP2O5 *= rnorm;
	nS *= rnorm;
	nCr2O3 *= rnorm;
	nZrO2 *= rnorm;
	nH2O *= rnorm;
	nCO2 *= rnorm;
	nCl *= rnorm;

	return 100.0;
}
