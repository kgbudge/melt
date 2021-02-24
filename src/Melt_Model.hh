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

#include "ds++/Assert.hh"

#include "State.hh"

//-----------------------------------------------------------------------------//
class Melt_Model : public State_Base
{
	public:

		struct Reaction
        {
          unsigned i; // extracted phase
          double dGf0; // initial rate of change of free energy of reaction
          double dGfs; // rate of change of free energy of solid phase with reaction
          double extent; // typically, the extent to which the reaction can proceed before hitting a constraint.

          unsigned nz;  // Number of phases in reaction
          unsigned p[E_END]; // Phases in reaction
          double n[E_END]; // Phase amount in reaction (negative == reagent)
        };

	    Melt_Model(double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
                   double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
                   double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
                   double nFe2O3,  double nZrO2);

	    explicit Melt_Model(double const Z[E_END]);

		double Z(unsigned i) const { Require(i<E_END); return Z_[i]; }

		double Gf(double const XP[P_END]) const;
		 
		double dGf(double const XP[P_END], 
 				   double const xm[E_END],
                   double Gf, 
                   unsigned phase) const;
		 
		double dGfm(double const XP[P_END], 
 				    double const xm[E_END],
                    double Gf, 
                    unsigned phase) const;

		double Gfm(double const XM[E_END]) const;

		double Gfmelt(double const X[P_END],
	                  unsigned NMP,
	                  Reaction const cphase[P_END], 
                      double p[M_END],
                      double e);
		
		Phase minimize_Gf(double XP[P_END]);

        void print(Reaction const &) const;
        void print_reverse(Reaction const &) const;

	private:

// **** Implementation

        double calculate_extent_(Reaction const &, double const XP[P_END], double const xm[E_END], double dir) const;
        void compute_current_melt_composition_(double const XP[P_END], double xm[E_END]) const;
        void compute_current_solid_composition_(double const XP[P_END], double xs[E_END]) const;
		void compute_melt_endmember_composition_(double const xm[E_END], double xend[M_END]) const;

	    Reaction construct_reaction_(double const XP[],
								     double const xm_[],
								     double const xs_[], 
									 double Gf, 
                                     double Gfm, 
									 unsigned i) const;

	    double minimize_trial_set_(unsigned NC,
							       Reaction const cphase[P_END],
								   double XP[P_END]);

        void update_current_state_(double XP[P_END],
							       double xm[E_END],
							       double xs[E_END]) const;

        void update_solid_state_(double XP[P_END]) const;

// **** Data

        // Fixed attributes of melt
		double Z_[E_END]; // amount of each fusible element
		double cnorm_; // total moles of atoms; used for normalization of various tests
};

#endif // melt_hh
