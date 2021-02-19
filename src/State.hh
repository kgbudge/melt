/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * State.hh
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

#ifndef State_HH
#define State_HH

#include <memory>
#include <string>
#include <vector>

#include "Phase.hh"
#include "element.hh"

// Describes a state of a sample of geological material. 
class State
{
  public:
    // Construct a starting state using the usual composition by oxide.
    // Each of these phases is in the standard phase table and the starting
    // state will have these as the active phases.
	State(std::string const &name,     
		  double T,
		  double P,   
          double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
          double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
          double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
          double nFe2O3,  double nZrO2);

    // Construct a starting state from moles of each bare element.
    // This will be translated into starting oxide phases.
	State(std::string const &name,
          double T,
          double P,
          double const x[E_END]);

   // Construct a state congruent with a Melt_Model and with the given starting phases
	State(std::string const &name,
		  class Melt_Model const &,
          unsigned NP,
		  unsigned const cphase[],
	      double const xphase[]);

// Accessors

	constexpr double T() const { return T_; }
	constexpr double P() const { return P_; }
	constexpr std::vector<Phase> const &phase() const { return phase_; }
	constexpr std::vector<double> const &Gf() const { return Gf_; }
	constexpr int const *ph() const { return ph_; }
	constexpr double const *X() const { return X_; }
    constexpr bool const *is_element_active() const { return is_element_active_; }
	constexpr double const *V() const { return V_; }

// Services

    // Determine the set of solid phases that minimizes the free energy of the 
    // sample.
    bool do_ladder_update();

    void update();

    Phase melt() const;

  private:

    void compute_gf_();

    std::string name_;         // Name of the state, usually something coded like "1.2.3"
    double T_, P_;             // ambient conditions
	std::vector<Phase> phase_; // Phases from which state is constructed
	std::vector<double> Gf_;   // Free energy of each phase
	int ph_[E_END];            // Phase index of each phase present in the sample.
	double X_[E_END];          // Number of moles of each phase in the sample.
    bool is_element_active_[E_END];  // Which chemical elements are present?
	double element_activity_[E_END]; // Current element activities
	double V_[E_END];          // Molar volume of each active phase

};

#endif // State_HH
