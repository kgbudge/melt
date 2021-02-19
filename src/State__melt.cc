/*
 * melt.cc
 * Copyright (C) 2019 Kent G. Budge <kgb@kgbudge.com>
 * 
 * norm is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * norm is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "State.hh"

#include "Melt_Model.hh"

//-----------------------------------------------------------------------------//
Phase State::melt() const
{
	using namespace std;
	
	Melt_Model model(*this);

    // Initial guess is that all fusible phases are fully melted. 
	// We will then see what should crystallize out.

	double XP[P_END];
	fill(XP, XP+P_END, 0.0);
	
	return model.minimize_Gf(XP);

}
