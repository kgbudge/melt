/*
 * Melt_Model__mimimize_Gf.cc
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

#include <iomanip>
#include <iostream>
#include <limits>
#include <set>

#include "Model.hh"

//-----------------------------------------------------------------------------//
/*! Find the mixture of solid phases and a melt phase that minimized total Gibbs
 * free energy.
 * 
 * \param[in,out] XP Contains initial guess of number of moles of each solid 
 * phase. This is compressed to include only the Melt::NP_ active phases 
 * specified by Melt::p_. On return, contains number of moles of each solid 
 * phase that minimizes free energy.
 *
 * \return The melt phase that minimizes free energy.
 */
Phase Melt_Model::minimize_Gf(double XP[P_END])
{
	using namespace std;

	double xm[E_END];  // Current estimate of melt composition

	unsigned const NPL = phase_.size();  // Number of solid phases in phase library.

	// Restart loop. This discards the search space and begins building up a
	// new search space for the minimum of free energy.
	for (;;)
	{
		unsigned cphase[E_END];
		unsigned NMP = 0;
		compute_current_melt_composition_(XP, xm);
		
		// Main loop. Here we build up a search space for the minimum of free energy.
		for (;;)
		{
			// Calculate free energy of current estimate of sample state, and its
			// gradient with every solid phase in library.

			double Gf = this->Gf(XP);
			cout << "Gf = " << fixed << setprecision(3) << Gf << endl;
			double g[P_END];
			bool has_crystallizing = false, has_melting = false;
			unsigned NS = 0; // count of number of solid phases
			for (unsigned i=0; i<NPL; ++i)
			{
				g[i] = dGf(XP, xm, Gf, i);
				if (g[i]<-1.0e-6*cnorm_)
				{
					has_crystallizing = true;
				}
				else if (g[i]>1.0e-6*cnorm_ && XP[i]>1.0e-9*cnorm_)
				{
					has_melting = true;
				}
				if (XP[i]>1.0e-9*cnorm_)
				{
					NS++;
				}
				else
				{
					XP[i] = 0; // Prune trace phase
				}
			}

			// If nothing can change, we must be done.

			if (has_crystallizing || has_melting)
			{
				if (has_melting)
				{
					// There are spontaneously melting phases. Do them first.
					NMP = 0;
					cout << "Melting phases:" << endl;
					for (unsigned i=0; i<NPL; ++i)
					{
						if (g[i]>1.0e-6*cnorm_ && XP[i]>1.0e-9*cnorm_) 
						{
							cphase[NMP++] = i;	
							cout << "  " << phase_[i].name << " g = " << g[i] << endl;
						}
					}
				}
				else
				{
					// Update solid phases
					double xs[E_END];
					compute_current_solid_composition_(XP, xs);
					set<unsigned> in_solid;
					{
						State solid_state("s", T_, P_, xs);
						solid_state.do_ladder_update();
						auto const &solid_state_ph = solid_state.ph();
						auto const &solid_state_X = solid_state.X();
						auto const &solid_state_phase = solid_state.phase();
						for (unsigned i=0; i<E_END; ++i)
						{
							if (solid_state_X[i]>0.0)
							{
								in_solid.insert(solid_state_phase[solid_state_ph[i]].index);
							}
						}
					}

					// Prune anything from the old list that is exhausted.
					cout << "Crystallization phases:" << endl;
					unsigned n = 0;
					set<unsigned> cset;
					for (unsigned i=0; i<NMP; ++i)
					{
						Phase const &phase = phase_[cphase[i]];
						unsigned const N = phase.nz;
						bool keep = true;
						for (unsigned j=0; j<N; ++j)
						{
							if (xm[phase.z[j]]<1.0e-9*cnorm_)
							{
								keep = false;
								break;
							}
						}
						if (keep)
						{
							cout << "  " << phase_[cphase[i]].name << " g = " << g[cphase[i]] << endl;
							cset.insert(cphase[i]);
							cphase[n++] = cphase[i];	
						}
					}
					NMP = n;


					// We have possible crystallizing phases. For each such phase,
					// see if the phase requires any element not in the mix. If so,
					// we will try to construct an incongruent reaction to produce it
					// and see if that is spontaneous.

					double merit = 0.0;
					for (unsigned i=0; i<NPL; ++i)
					{
						if (cset.count(i) == 0 && g[i]<0.0) 
						{
							Phase const &ph = phase_[i];
							unsigned N = ph.nz; 
							double xsn[E_END];
							copy(xs, xs+E_END, xsn);
							double x1 = std::numeric_limits<double>::max();
							for (unsigned j=0; j<N; ++j)
							{
								unsigned z = ph.z[j]; 
								x1 = min(x1, xm[z]/ph.n[j]);
								xsn[z] += ph.n[j]*0.001*cnorm_;
							}
							if (x1>1.0e-9*cnorm_)
							{
								State state("t",  T_, P_, xsn);
								state.do_ladder_update();
								bool found = false;
								auto const &state_ph = state.ph();
								auto const &state_phase = state.phase();
								auto const &state_x = state.X();
								for (unsigned j=0; j<E_END; j++)
								{
									if (state_phase[state_ph[j]].index==phase_[i].index)
									{
										found = state_x[j]>0.0; 
										break;
									}
								}
								if (found)
								{
									double mm = -g[i]*x1;
									cout << "  " << phase_[i].name << " g = " << mm << endl;
									if  (mm>merit)
									{
										merit =  mm;
										cphase[NMP] = i;
									}
								}
							}
							else
							{
								// Can only crystallize incongruously
								//						Insist(false, "construction");
							}
						}
					}
					NMP++;
				}

				// Add a single phase that is our best guess of the next melt element

				// Minimize on this melt set.

				double revised_Gf = minimize_trial_set_(NMP, cphase, XP);
				compute_current_melt_composition_(XP, xm);

				double mtot = 0.0;
				for (unsigned i=0; i<E_END; ++i)
				{
					if (xm[i]<1.0e-9*cnorm_)
					{
						xm[i] = 0.0;
					}
					mtot += xm[i];
				}	
				if (mtot<=0.0 || Gf - revised_Gf  < 1.0e-9*cnorm_)
				{
					break;
				}
			}
		}
	}
		
	double V = 0.0;
	double xend[M_END];
	compute_melt_endmember_composition_(xm, xend);
	for (unsigned i=0; i<M_END; ++i)
	{
		double xi = xend[i];
		if (i != E_O && xi>0.0)
		{
			Phase const &phase = ::phase[melt_endmember[i]];
			V += xi*phase.model->volume(phase, T_, P_);
		}
	}
		
	Phase Result;
	Result.index = 0;
	Result.name = "melt";
	Result.nz = 0;
	Result.V = V;
	for (unsigned i=0; i<E_END; ++i)
	{
		if (xm[i]>1e-9)
		{
			Result.z[Result.nz] = i;
			Result.n[Result.nz] = xm[i];
			Result.nz++;
		}
	}
	// If Result.nz is zero, this flags to client that there is no melt.

	Result.S0 = 0.0;
	Result.Hf0 = Gfm(xm);
	Result.model = MELT;

	return Result;
}
