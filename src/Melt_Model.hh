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

enum Melt_Endmember
{
		M_H2O,
		M_SiO2,
		M_Al2O3,
		M_Mg2SiO4,
		M_Fe2SiO4,
		M_CaSiO3,
		M_Na2SiO3,
		M_KAlSi2O6,

		ME_END,

        M_CO2 = ME_END,  // 8
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

unsigned const melt_endmember[M_END] =
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

//-----------------------------------------------------------------------------//
class Melt_Model
{
	public:
		explicit Melt_Model(State const &state) noexcept(false);

		unsigned NP() const noexcept { return NP_; }
		double Z(unsigned i) const { Require(i<E_END); return Z_[i]; }

		double Gf(double const XP[P_END]) const;
		 
		double dGf(double XP[P_END], double Gf0, unsigned phase) const;

		double Gfm(double const XM[E_END]) const;

		double Gfmelt(double const X[P_END], double p[M_END], double e) const;
		
		Phase minimize_Gf(double XP[P_END]);

	private:

        void compute_current_melt_composition_(double const XP[], double xm[]) const;
	    double minimize_trial_set_(double XP[P_END]);
			 
        // Copied from parent State 
        double T_;
	    double P_;
	    std::vector<Phase> phase_; // Phases from which state is constructed
	    std::vector<double> Gf_;   // Free energy of each phase

        // Fixed attributes of melt
	    std::vector<bool> is_fusible_;  // Is this phase fusible?
		double Z_[E_END]; // amount of each fusible element
		double Gfr_; // non-meltable phases total free energy
		double Gfm_[M_END]; // Melt end member free energies

        // Current minimization set. Must include all fusible solid phases in
        // the trial state plus the trial set of crystallizing phases.
		unsigned NP_; // number of active phases
		unsigned ph_[P_END]; // indices of active phases
};

#endif // melt_hh
