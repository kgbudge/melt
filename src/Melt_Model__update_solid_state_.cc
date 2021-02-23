/*
 * Melt_Model__update_solid_state_.cc
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
/*! Do a ladder update on the solid phases. This converts them to the minimum
 * free energy configuration ignoring melt.
 * 
 * \param[in]  XP Solid phase composition
  */
void Melt_Model::update_solid_state_(double XP[P_END])
{
	using namespace std;
	
	compute_current_solid_composition_(XP);
	State solid("s", T_, P_, xs_);
	solid.do_ladder_update();
	fill(XP, XP+NP_, 0.0);
	auto const &solid_X = solid.X();
	auto const &solid_phase = solid.phase();
	auto const &solid_ph = solid.ph();
	for (unsigned i=0; i<E_END; ++i)
	{
		if (solid_X[i]>1.0e-6*cnorm_)
		{
			unsigned gi = solid_phase[solid_ph[i]].index;
			Check(gi<P_END);
			gi = imap_[gi];
			Check(gi<NP_);
			XP[gi] = solid_X[i];
		}
	}
}
