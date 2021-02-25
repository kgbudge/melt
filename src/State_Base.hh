/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * State_Base.hh
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

#ifndef State_Base_HH
#define State_Base_HH

#include "gsl/gsl_linalg.h"

#include "Phase.hh"
#include "element.hh"
#include "melt_model.hh"
#include "phase_enum.hh"

// Shared static libraries and data structures for State and Melt_Model
class State_Base
{
  public:
 
// Static

    static 
    void initialize_globals(double T, double P,   
          double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
          double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
          double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
          double nFe2O3,  double nZrO2);

	static double T() { return T_; }
	static double P() { return P_; }
	static unsigned NP() { return NP_; }
	static constexpr Phase const *phase() { return phase_; }
	static constexpr double const *Gf() { return Gf_; }

  protected:

// Static
    static double T_, P_;                   // ambient conditions
    static bool is_element_active_[E_END];  // Which chemical elements are present?
	static unsigned NP_;                    // Number of phases in phase library
	static Phase phase_[P_END];             // Phase library
	static double Gf_[P_END];               // Free energy of each phase in phase library
	static bool is_fusible_[P_END];         // Is this phase fusible?
    static double Gfm_[M_END];              // Melt endmember free energies

	static gsl_matrix *gsl_A_;
	static gsl_matrix *gsl_V_;
	static gsl_vector *gsl_S_;
	static gsl_vector *gsl_work_;
	static gsl_vector *gsl_b_;
	static gsl_vector *gsl_x_;	
	static gsl_vector *gsl_aGf_;
	static gsl_permutation *gsl_permute_;
};

#endif // State_HH
