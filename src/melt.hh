// melt.hh
//
// Copyright (C) 2019 - Kent G. Budge
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef melt_hh
#define melt_hh

#include <vector>

#include "element.hh"
#include "phase_enum.hh"
#include "State.hh"

enum Endmember
{
		M_H2O,
		M_SiO2,
		M_Al2O3,
		M_Mg2SiO4,
		M_Fe2SiO4,
		M_CaSiO3,
		M_Na2SiO3,
		M_KAlSi2O6,

        M_CO2,  // 8
        M_NaAlSiO4,
        M_NaAlSi3O8,
        M_MgO,
        M_S2,
	    M_NaCl,
        M_Mg2Si2O6, //14
		M_CaO,
        M_CaMgSi2O6,
        M_CaAl2Si2O8, // 17
	    M_KAlSi3O8,

        M_END
};

unsigned const endmember_element[M_END] =
{
		E_H,
		E_SI,
		E_AL,
		E_MG,
		E_FE,
		E_CA,
		E_NA,
		E_K
};

unsigned const endmember[M_END] =
{
		P_WATER_VAPOR,
		P_SiO2_LIQUID,
		P_CORUNDUM_LIQUID,
		P_FORSTERITE_LIQUID,
		P_FAYALITE_LIQUID,
		P_WOLLASTONITE_LIQUID,
		P_Na2SiO3_LIQUID,
		P_LEUCITE_LIQUID,

        P_CO2,
        P_NEPHELINE_LIQUID,
        P_ALBITE_LIQUID,
        P_PERICLASE_LIQUID,
        P_DISULFUR,
	    P_HALITE_LIQUID,
        P_ENSTATITE_LIQUID,
		P_LIME_LIQUID,
        P_DIOPSIDE_LIQUID,
        P_ANORTHITE_LIQUID,
	    P_K_FELDSPAR_LIQUID,
};

double const mixN[M_END] = 
{
		1, // M_H2O: We expect water to really reduce the entropy
		0.25, // M_SiO2,
		0.5, // M_Al2O3,
		1, // M_Mg2SiO4,
		1, // M_Fe2SiO4,
		0.5, // M_CaSiO3,
		2, // M_Na2SiO3,
		1, // M_KAlSi2O6,

        1, // M_CO2,
        1, // M_NaAlSiO4,
        1, // M_NaAlSi3O8,
        1, // M_MgO,
        1, // M_S2,
	    1, // M_NaCl,
        1, // M_Mg2Si2O6,
		1, // M_CaO,
        1, // M_CaMgSi2O6,
        1, // M_CaAl2Si2O8,
	    1, // M_KAlSi3O8,
};

double melt(double T, 
          double P, 
          std::vector<struct Phase> const &phase,
	      std::vector<double> const &Gf,
          State const &state, 
          struct Phase &new_phase);

#endif // melt_hh
