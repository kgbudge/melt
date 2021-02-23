/*
 * Melt_Model__update_current_state.cc
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

//-----------------------------------------------------------------------------//
/*! Update the current state from the specified solid phase composition.
 *
 * \param XP Current solid phase composition.
 */ 
void Melt_Model::update_current_state_(double XP[P_END],
                                       double xm[E_END],
                                       double xs[E_END]) const
{
	// Do a ladder update to minimize the solid phases
	update_solid_state_(XP);
	compute_current_melt_composition_(XP, xm);
	compute_current_solid_composition_(XP, xs);
};

