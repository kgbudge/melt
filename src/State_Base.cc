/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * State_Base.cc
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

#include "State_Base.hh"

// static
double State_Base::T_;
double State_Base::P_;
bool State_Base::is_element_active_[E_END];  // Which chemical elements are present?
unsigned State_Base::NP_;                    // Number of phases in phase library
Phase State_Base::phase_[P_END];             // Phase library
double State_Base::Gf_[P_END];               // Free energy of each phase in phase library
bool State_Base::is_fusible_[P_END];         // Is this phase fusible?
double State_Base::Gfm_[M_END];              // Melt endmember free energies

gsl_matrix *State_Base::gsl_A_;
gsl_matrix *State_Base::gsl_V_;
gsl_vector *State_Base::gsl_S_;
gsl_vector *State_Base::gsl_work_;
gsl_vector *State_Base::gsl_b_;
gsl_vector *State_Base::gsl_x_;	
gsl_vector *State_Base::gsl_aGf_;
gsl_permutation *State_Base::gsl_permute_;
