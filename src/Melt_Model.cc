/*
 * Melt_Model.cc
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

#include "Melt_Model.hh"

#include "Model.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

//-----------------------------------------------------------------------------//
/*! Create a Melt_Model reflecting possible melting of an existing State.
 *
 * \param state Current state of the mineral ensemble.
 */ 
Melt_Model::Melt_Model(double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
                       double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
                       double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
                       double nFe2O3,  double nZrO2)
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

	// Compute the fully melted composition (excluding phases for which we
	// do not have a liquid counterpart). Determine free energy contribution
	// of these nonfusible phases.

	Z_[E_H] = 2*nH2O;
	Z_[E_C] = nCO2;
	Z_[E_NA] = 2*nNa2O;
	Z_[E_MG] = nMgO;
	Z_[E_AL] = 2*nAl2O3;
	Z_[E_SI] = nSiO2;
	Z_[E_P] = 2*nP2O5;
	Z_[E_S] = nS;
	Z_[E_CL] = nCl;
	Z_[E_K] = 2*nK2O;
	Z_[E_CA] = nCaO;
	Z_[E_TI] = nTiO2;
	Z_[E_CR] = 2*nCr2O3;
	Z_[E_MN] = nMnO;
	Z_[E_FE] = nFeO + 2*nFe2O3;
	Z_[E_ZR] = nZrO2;
	Z_[E_O] = nH2O + 2*nCO2 + nNa2O + nMgO + 3*nAl2O3 + 2*nSiO2 + 5*nP2O5 +
		2*nK2O + nCaO + 2*nTiO2 + 3*nCr2O3 + nMnO + nFeO + 3*nFe2O3 + 2*nZrO2;

	cout << "Full melt elemental molar composition:" << endl << defaultfloat;
	cnorm_ = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		Check(Z_[i]>=0.0);
		cnorm_ += Z_[i];
		if (Z_[i]>0)
		{
			cout << element_name[i] << ": " << setprecision(3) << Z_[i] << endl;
		}
	}
}


//-----------------------------------------------------------------------------//
Melt_Model::Melt_Model(double const Z[E_END])
{
	using namespace std;

	for (unsigned i=0; i<E_END; ++i)
	{
		Require(Z[i]>=0.0);
		Z_[i] = Z[i];
	}

	cout << "Full melt elemental molar composition:" << endl << defaultfloat;
	cnorm_ = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		Check(Z_[i]>=0.0);
		cnorm_ += Z_[i];
		if (Z_[i]>0)
		{
			cout << element_name[i] << ": " << setprecision(3) << Z_[i] << endl;
		}
	}
}

char const * const endmember_element_name[] =
{
	"H",
	"Si",
	"Al",
	"Mg",
	"Fe(+2)",
	"Ca",
	"Na",
	"K"
};

