/*
 * Melt_Model__construct_reaction_.cc
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
/*! Characterize the consequences of moving a specified solid phaes composition
 * from the melt to the solid phase.
 *
 * \param XP Current solid phase composition
 */
 Melt_Model::Reaction Melt_Model::construct_reaction_(double const XP[],
                                                      double const xm[],
                                                      double const xs[], 
                                                      double const Gf, 
                                                      double const Gfm, 
                                                      unsigned const i) const
{
	 using namespace std;

	 Reaction Result;
	 Result.i = i;

	 // Compute rate of change of free energy of melt phase
	 double const dGfm = this->dGfm(XP, xm, Gfm, i);

	 // Compute rate of change of free energy of solid phase

	 Phase const &ph = phase_[i];
	 unsigned const N = ph.nz;
	 double xsn[E_END];
	 copy(xs, xs+E_END, xsn);
	 for (unsigned j=0; j<N; ++j)
	 {
		 unsigned z = ph.z[j]; 
		 xsn[z] += ph.n[j]*0.001*cnorm_;
	 }
	 State state(ph.name, xsn);
	 state.do_ladder_update();
	 auto const &state_ph = state.ph();
	 auto const &state_phase = state.phase();
	 auto const &state_x = state.X();
	 double XPP[P_END];
	 fill(XPP, XPP+P_END, 0.0);
	 for (unsigned j=0; j<E_END; ++j)
	 {
		 if (state_x[j]>0.0)
		 {
			 XPP[state_ph[j]] = state_x[j];
		 }
	 }
	 double dGfs = 0.0;
	 Result.nz = 0;
	 for (unsigned i=0; i<NP_; ++i)
	 {
		 if (fabs(XP[i] - XPP[i])>1.0e-12*cnorm_)
		 {
			 double dn = (XPP[i] - XP[i])*1000./cnorm_;
			 dGfs += dn*Gf_[i];
			 Result.n[Result.nz] = dn;
			 Result.p[Result.nz] = i;
			 Result.nz++;
		 }
	 }
	 Result.dGfs = dGfs;
	 Result.dGf0 = dGfs + dGfm;

	 Result.extent = calculate_extent_(Result, XP, xm, Result.dGf0);
	 return Result;
}