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

#include <numeric>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"

#include "ds++/Assert.hh"

#include "algorithm.hh"
#include "constants.hh"
#include "D1.hh"
#include "element.hh"
#include "Melt_Model.hh"
#include "Model.hh"
#include "Phase.hh"
#include "State.hh"

//-----------------------------------------------------------------------------//
using namespace std;

//-----------------------------------------------------------------------------//
double const Wnorm = 1.0e-3;
double const Ww_q = 0.*Wnorm;
double const Ww_cor = -0.*Wnorm;
double const Ww_fo = 0.*Wnorm;
double const Ww_fa = 0.*Wnorm;
double const Ww_wo = 0.*Wnorm;
double const Ww_sm = -0.*Wnorm;
double const Ww_kal = 0.*Wnorm;
double const Ww_co2 = 0*Wnorm;
double const Ww_ne = 0*Wnorm;
double const Ww_ab = 0*Wnorm;
double const Ww_pc = 0*Wnorm;
double const Ww_s = 0*Wnorm;
double const Ww_ha = 0*Wnorm;
double const Ww_en = 0*Wnorm;
double const Ww_li = 0*Wnorm;
double const Ww_mdi = 0*Wnorm;
double const Ww_an = 0*Wnorm;
double const Ww_or = 0*Wnorm;
double const Wq_cor = -1230*Wnorm; // calibrated from SiO2-Al2O3 eutectic at 1868K=1595C  Al 0.9m vs. 0.4m
double const Wq_fo = 0; // do not coexist
double const Wq_fa = 2400.*Wnorm; // calibrated from Q-Fa eutectic near 1422K=1149C
double const Wq_wo = 1850.*Wnorm;// calibrated from eutectic of 1699K=1426C
double const Wq_sm = -3050*Wnorm;// calibrated from eutectic at 1062K=789C
double const Wq_kal = -0.*Wnorm;
double const Wq_co2 = 0*Wnorm;
double const Wq_ne = 0*Wnorm;
double const Wq_ab = 0*Wnorm;
double const Wq_pc = 0*Wnorm;
double const Wq_s = 0*Wnorm;
double const Wq_ha = 0*Wnorm;
double const Wq_en = 0*Wnorm;
double const Wq_li = 0*Wnorm;
double const Wq_mdi = 0*Wnorm;
double const Wq_an = 0*Wnorm;
double const Wq_or = 0*Wnorm;
double const Wcor_fo = 0*Wnorm;
double const Wcor_fa = 0*Wnorm;
double const Wcor_wo = 0*Wnorm;
double const Wcor_sm = 0*Wnorm;
double const Wcor_kal = 0*Wnorm;
double const Wcor_co2 = 0*Wnorm;
double const Wcor_ne = 0*Wnorm;
double const Wcor_ab = 0*Wnorm;
double const Wcor_pc = 0*Wnorm;
double const Wcor_s = 0*Wnorm;
double const Wcor_ha = 0*Wnorm;
double const Wcor_en = 0*Wnorm;
double const Wcor_li = 0*Wnorm;
double const Wcor_mdi = 0*Wnorm;
double const Wcor_an = 0*Wnorm;
double const Wcor_or = 0*Wnorm;
double const Wfo_fa = -0.*Wnorm;
double const Wfo_wo = -0.*Wnorm;
double const Wfo_sm = -0.*Wnorm;
double const Wfo_kal = 0.*Wnorm;
double const Wfo_co2 = 0*Wnorm;
double const Wfo_ne = 0*Wnorm;
double const Wfo_ab = -6000*Wnorm; // calibrated from eutectic at 1376K=1103C
double const Wfo_pc = 0*Wnorm;
double const Wfo_s = 0*Wnorm;
double const Wfo_ha = 0*Wnorm;
double const Wfo_en = 9000*Wnorm; // calibrated from Q-Fo eutectic at 1815K=1542C
double const Wfo_li = 0*Wnorm;
double const Wfo_mdi = 30000*Wnorm; // calibrated from eutectic at 1384C
double const Wfo_an = -32500*Wnorm;     // calibrated from eutectic at 1317C (with spinel formation)
double const Wfo_or = 0*Wnorm;
double const Wfa_wo = -0.*Wnorm;
double const Wfa_sm = -0.*Wnorm;
double const Wfa_kal = 0.*Wnorm;
double const Wfa_co2 = 0*Wnorm;
double const Wfa_ne = 0*Wnorm;
double const Wfa_ab = 0*Wnorm;
double const Wfa_pc = 0*Wnorm;
double const Wfa_s = 0*Wnorm;
double const Wfa_ha = 0*Wnorm;
double const Wfa_en = 0*Wnorm;
double const Wfa_li = 0*Wnorm;
double const Wfa_mdi = 0*Wnorm;
double const Wfa_an = 0*Wnorm;
double const Wfa_or = 0*Wnorm;
double const Wwo_sm = -0.*Wnorm;
double const Wwo_kal = 0.*Wnorm;
double const Wwo_co2 = 0*Wnorm;
double const Wwo_ne = 0*Wnorm;
double const Wwo_ab = 0*Wnorm;
double const Wwo_pc = 0*Wnorm;
double const Wwo_s = 0*Wnorm;
double const Wwo_ha = 0*Wnorm;
double const Wwo_en = 0*Wnorm;
double const Wwo_li = 0*Wnorm;
double const Wwo_mdi = 0*Wnorm;
double const Wwo_an = 0*Wnorm;
double const Wwo_or = 0*Wnorm;
double const Wsm_kal = 0.*Wnorm;
double const Wsm_co2 = 0*Wnorm;
double const Wsm_ne = 0*Wnorm;
double const Wsm_ab = 0*Wnorm;
double const Wsm_pc = 0*Wnorm;
double const Wsm_s = 0*Wnorm;
double const Wsm_ha = 0*Wnorm;
double const Wsm_en = 0*Wnorm;
double const Wsm_li = 0*Wnorm;
double const Wsm_mdi = 0*Wnorm;
double const Wsm_an = 0*Wnorm;
double const Wsm_or = 0*Wnorm;
double const Wkal_co2 = 0*Wnorm;
double const Wkal_ne = 0*Wnorm;
double const Wkal_ab = 0*Wnorm;
double const Wkal_pc = 0*Wnorm;
double const Wkal_s = 0*Wnorm;
double const Wkal_ha = 0*Wnorm;
double const Wkal_en = 0*Wnorm;
double const Wkal_li = 0*Wnorm;
double const Wkal_mdi = 0*Wnorm;
double const Wkal_an = 0*Wnorm;
double const Wkal_or = 0*Wnorm;
double const Wco2_ne = 0*Wnorm;
double const Wco2_ab = 0*Wnorm;
double const Wco2_pc = 0*Wnorm;
double const Wco2_s = 0*Wnorm;
double const Wco2_ha = 0*Wnorm;
double const Wco2_en = 0*Wnorm;
double const Wco2_li = 0*Wnorm;
double const Wco2_mdi = 0*Wnorm;
double const Wco2_an = 0*Wnorm;
double const Wco2_or = 0*Wnorm;
double const Wne_ab = 0*Wnorm;
double const Wne_pc = 0*Wnorm;
double const Wne_s = 0*Wnorm;
double const Wne_ha = 0*Wnorm;
double const Wne_en = 0*Wnorm;
double const Wne_li = 0*Wnorm;
double const Wne_mdi = 0*Wnorm;
double const Wne_an = 0*Wnorm;
double const Wne_or = 0*Wnorm;
double const Wab_pc = 0*Wnorm;
double const Wab_s = 0*Wnorm;
double const Wab_ha = 0*Wnorm;
double const Wab_en = 0*Wnorm;
double const Wab_li = 0*Wnorm;
double const Wab_mdi = 0*Wnorm;
double const Wab_an = 0*Wnorm;
double const Wab_or = 0*Wnorm;
double const Wpc_s = 0*Wnorm;
double const Wpc_ha = 0*Wnorm;
double const Wpc_en = 0*Wnorm;
double const Wpc_li = 0*Wnorm;
double const Wpc_mdi = 0*Wnorm;
double const Wpc_an = 0*Wnorm;
double const Wpc_or = 0*Wnorm;
double const Ws_ha = 0*Wnorm;
double const Ws_en = 0*Wnorm;
double const Ws_li = 0*Wnorm;
double const Ws_mdi = 0*Wnorm;
double const Ws_an = 0*Wnorm;
double const Ws_or = 0*Wnorm;
double const Wha_en = 0*Wnorm;
double const Wha_li = 0*Wnorm;
double const Wha_mdi = 0*Wnorm;
double const Wha_an = 0*Wnorm;
double const Wha_or = 0*Wnorm;
double const Wen_li = 0*Wnorm;
double const Wen_mdi = 0*Wnorm;
double const Wen_an = 0*Wnorm; 
double const Wen_or = 0*Wnorm;
double const Wli_mdi = 0*Wnorm;
double const Wli_an = 0*Wnorm;
double const Wli_or = 0*Wnorm;
double const Wmdi_an = -13000*Wnorm; // Calibrated from eutectic at 1274C
double const Wmdi_or = 0*Wnorm;

double const W[M_END][M_END] =
{
	{0.0,      Ww_q,  Ww_cor,    Ww_fo,    Ww_fa,    Ww_wo,   Ww_sm,   Ww_kal,   Ww_co2,   Ww_ne,     Ww_ab,   Ww_pc,   Ww_s,   Ww_ha,   Ww_en,   Ww_li,   Ww_mdi,   Ww_an, Ww_or},
	{Ww_q,      0.0,  Wq_cor,    Wq_fo,    Wq_fa,    Wq_wo,   Wq_sm,   Wq_kal,   Wq_co2,   Wq_ne,     Wq_ab,   Wq_pc,   Wq_s,   Wq_ha,   Wq_en,   Wq_li,   Wq_mdi,   Wq_an, Wq_or},
	{Ww_cor, Wq_cor,     0.0,  Wcor_fo,  Wcor_fa,  Wcor_wo, Wcor_sm, Wcor_kal, Wcor_co2, Wcor_ne,   Wcor_ab, Wcor_pc, Wcor_s, Wcor_ha, Wcor_en, Wcor_li, Wcor_mdi, Wcor_an, Wcor_or},
    {Ww_fo,   Wq_fo, Wcor_fo,      0.0,   Wfo_fa,   Wfo_wo,  Wfo_sm,  Wfo_kal,  Wfo_co2,  Wfo_ne,    Wfo_ab,  Wfo_pc,  Wfo_s,  Wfo_ha,  Wfo_en,  Wfo_li,  Wfo_mdi,  Wfo_an, Wfo_or},
    {Ww_fa,   Wq_fa, Wcor_fa,   Wfo_fa,      0.0,   Wfa_wo,  Wfa_sm,  Wfa_kal,  Wfa_co2,  Wfa_ne,    Wfa_ab,  Wfa_pc,  Wfa_s,  Wfa_ha,  Wfa_en,  Wfa_li,  Wfa_mdi,  Wfa_an, Wfa_or},
    {Ww_wo,   Wq_wo, Wcor_wo,   Wfo_wo,   Wfa_wo,      0.0,  Wwo_sm,  Wwo_kal,  Wwo_co2,  Wwo_ne,    Wwo_ab,  Wwo_pc,  Wwo_s,  Wwo_ha,  Wwo_en,  Wwo_li,  Wwo_mdi,  Wwo_an, Wwo_or},
    {Ww_sm,   Wq_sm, Wcor_sm,   Wfo_sm,   Wfa_sm,   Wwo_sm,     0.0,  Wsm_kal,  Wsm_co2,  Wsm_ne,    Wsm_ab,  Wsm_pc,  Wsm_s,  Wsm_ha,  Wsm_en,  Wsm_li,  Wsm_mdi,  Wsm_an, Wsm_or},
    {Ww_kal, Wq_kal, Wcor_kal, Wfo_kal,  Wfa_kal,  Wwo_kal, Wsm_kal,      0.0, Wkal_co2, Wkal_ne,   Wkal_ab, Wkal_pc, Wkal_s, Wkal_ha, Wkal_en, Wkal_li, Wkal_mdi, Wkal_an, Wkal_or},
    {Ww_co2, Wq_co2, Wcor_co2, Wfo_co2,  Wfa_co2,  Wwo_co2, Wsm_co2, Wkal_co2,      0.0, Wco2_ne,   Wco2_ab, Wco2_pc, Wco2_s, Wco2_ha, Wco2_en, Wco2_li, Wco2_mdi, Wco2_an, Wco2_or},
    {Ww_ne,   Wq_ne, Wcor_ne,   Wfo_ne,   Wfa_ne,   Wwo_ne,  Wsm_ne,  Wkal_ne,  Wco2_ne,     0.0,    Wne_ab,  Wne_pc,  Wne_s,  Wne_ha,  Wne_en,  Wne_li,  Wne_mdi,  Wne_an, Wne_or},
    {Ww_ab,   Wq_ab, Wcor_ab,   Wfo_ab,   Wfa_ab,   Wwo_ab,  Wsm_ab,  Wkal_ab,  Wco2_ab,  Wne_ab,       0.0,  Wab_pc,  Wab_s,  Wab_ha,  Wab_en,  Wab_li,  Wab_mdi,  Wab_an, Wab_or},
    {Ww_pc,   Wq_pc, Wcor_pc,   Wfo_pc,   Wfa_pc,   Wwo_pc,  Wsm_pc,  Wkal_pc,  Wco2_pc,  Wne_pc,   Wab_pc,      0.0,  Wpc_s,  Wpc_ha,  Wpc_en,  Wpc_li,  Wpc_mdi,  Wpc_an, Wpc_or},
    {Ww_s,     Wq_s,  Wcor_s,    Wfo_s,    Wfa_s,    Wwo_s,   Wsm_s,   Wkal_s,   Wco2_s,   Wne_s,    Wab_s,    Wpc_s,    0.0,   Ws_ha,   Ws_en,   Ws_li,   Ws_mdi,   Ws_an, Ws_or},
    {Ww_ha,   Wq_ha, Wcor_ha,   Wfo_ha,   Wfa_ha,   Wwo_ha,  Wsm_ha,  Wkal_ha,  Wco2_ha,  Wne_ha,   Wab_ha,   Wpc_ha,  Ws_ha,     0.0,  Wha_en,  Wha_li,  Wha_mdi,  Wha_an, Wha_or},
    {Ww_en,   Wq_en, Wcor_en,   Wfo_en,   Wfa_en,   Wwo_en,  Wsm_en,  Wkal_en,  Wco2_en,  Wne_en,   Wab_en,   Wpc_en,  Ws_en,  Wha_en,     0.0,  Wen_li,  Wen_mdi,  Wen_an, Wen_or},
    {Ww_li,   Wq_li, Wcor_li,   Wfo_li,   Wfa_li,   Wwo_li,  Wsm_li,  Wkal_li,  Wco2_li,  Wne_li,   Wab_li,   Wpc_li,  Ws_li,  Wha_li,  Wen_li,     0.0,  Wli_mdi,  Wli_an, Wli_or},
    {Ww_mdi,  Wq_mdi, Wcor_mdi, Wfo_mdi, Wfa_mdi,  Wwo_mdi, Wsm_mdi, Wkal_mdi, Wco2_mdi, Wne_mdi,  Wab_mdi,  Wpc_mdi, Ws_mdi, Wha_mdi, Wen_mdi, Wli_mdi,      0.0, Wmdi_an, Wmdi_or},
    {Ww_an,   Wq_an,  Wcor_an,   Wfo_an,  Wfa_an,   Wwo_an,  Wsm_an,  Wkal_an,  Wco2_an,  Wne_an,   Wab_an,   Wpc_an,  Ws_an,  Wha_an,  Wen_an,  Wli_an,  Wmdi_an,     0.0, Wmdi_or},
    {Ww_or,   Wq_or,  Wcor_or,  Wfo_or,   Wfa_or,   Wwo_or,  Wsm_or,  Wkal_or,  Wco2_or,  Wne_or,   Wab_or,   Wpc_or,  Ws_or,  Wha_or,  Wen_or,  Wli_or,  Wmdi_or,  Wmdi_or,    0.0}
};

template<typename Real>
Real Nlog(Real const &x)
{
	if (x>0.0)
	{
		return x*log(x);
	}
	else
	{
		return to_Real<Real>(0.0, size(x));
	}
}

//-----------------------------------------------------------------------------//
/*! Create a Melt_Model reflecting possible melting of an existing State.
 *
 * \param state Current state of the mineral ensemble.
 */ 
 Melt_Model::Melt_Model(State const &state)
:
 T_(state.T()), P_(state.P()), phase_(state.phase()), Gf_(state.Gf()), 
 is_fusible_(state.Gf().size())
{
	// Compute the fully melted composition (excluding phases for which we
	// do not have a liquid counterpart). Determine free energy contribution
	// of these nonfusible phases.
	 
	Gfr_ = 0.0;  // Gibbs free energy contribution of nonfusible phases 
    fill(Z_, Z_+E_END, 0.0); // Prepare to accumulate sample composition.
	                         //This will also be the starting melt. 
	 
	cout << "Sample phases:" << endl;
	double const *const state_X = state.X();
	int const *const state_ph = state.ph();
	unsigned const NPL = phase_.size();
	for (unsigned i=0; i<E_END; ++i)
	{
		double const x = state_X[i];
		double xpr[E_END];  // to store composition of phase
		fill(xpr, xpr+E_END, 0.0); 
		if (x>0.0)  // Is this phase actually present?
		{
			unsigned p = state_ph[i];
			Phase const &ph = phase_[p];	
			cout << ph.name << endl;
			unsigned const N = ph.nz; // Number of elements in the phase
			double xO = 0.0;          // To accumulate oxygen balance of the phase
			bool fusible = true;
			for (unsigned j=0; j<N && fusible; ++j) 
			{
				double const xj = x*ph.n[j];  // Number of moles of element j in the phase.
				switch(ph.z[j])             // Switch on the element atomic number
				{
					case E_H:
						xO -= 0.5*xj;
						xpr[E_H] += xj;
						break; 

				    case E_C:
						xO -= 2*xj;
						xpr[E_C] += xj;
						break;
						
					case E_O:
						xO += xj;
						break;

					case E_NA:
						xO -= 0.5*xj;
						xpr[E_NA] += xj;
						break;

					case E_MG:
						xO -= xj;
						xpr[E_MG] += xj;
						break;

					case E_AL:
						xpr[E_AL] += xj;
						xO -= 1.5*xj;
						break;

					case E_SI:
						xpr[E_SI] += xj;
						xO -= 2*xj;
						break;

					case E_S:
						xpr[E_S] += xj;
						break;

					case E_CL:
						xpr[E_CL] += xj;
						xO += 0.5*xj;
						break;

					case E_K:
						xpr[E_K] += xj;
						xO -= 0.5*xj;
						break;

					case E_CA:
						xO -= xj;
						xpr[E_CA] += xj;
						break;

					case E_FE:
						xpr[E_FE] += xj;
						break;

					default:
						// non-fusible mineral
						cout << "phase " << ph.name << " cannot melt." << endl;
						fusible = false;
						break;
				}
			}
			if (!fusible || fabs(xO-xpr[E_FE])>1e-9) // at present, cannot handle ferric or oxidized sulfur melts
			{
				cout << "  Not fusible" << endl;
				Gfr_ += x*Gf_[p];
			}
			else
			{
				for (unsigned m=0; m<E_END; ++m)
				{
					Z_[m] += xpr[m];
				}
			}
		}
	}
	cout << "Full melt elemental molar composition:" << endl << defaultfloat;
	for (unsigned i=0; i<E_END; ++i)
	{
		if (Z_[i]>0)
	  	  cout << element_name[i] << ": " << Z_[i] << endl;
	}

	// Find all potentiall fusible phases
	 for (unsigned p=0; p<NPL; ++p)
	 {
		 Phase const &ph = phase_[p];	
		 cout << ph.name << endl;
		 unsigned const N = ph.nz; // Number of elements in the phase
		 double xO = 0.0;          // To accumulate oxygen balance of the phase
		 double xFe = 0.0;
		 bool fusible = true;
		 for (unsigned j=0; j<N && fusible; ++j) 
		 {
			 double const xj = ph.n[j];  // Number of moles of element j in the phase.
			 switch(ph.z[j])             // Switch on the element atomic number
			 {
				 case E_H:
					 xO -= 0.5*xj;
					 break; 

				 case E_C:
					 xO -= 2*xj;
					 break;

				 case E_O:
					 xO += xj;
					 break;

				 case E_NA:
					 xO -= 0.5*xj;
					 break;

				 case E_MG:
					 xO -= xj;
					 break;

				 case E_AL:
					 xO -= 1.5*xj;
					 break;

				 case E_SI:
					 xO -= 2*xj;
					 break;

				 case E_S:
					 break;

				 case E_CL:
					 xO += 0.5*xj;
					 break;

				 case E_K:
					 xO -= 0.5*xj;
					 break;

				 case E_CA:
					 xO -= xj;
					 break;

				 case E_FE:
					 xFe += xj;
					 break;

				 default:
					 // non-fusible mineral
					 cout << "phase " << ph.name << " cannot melt." << endl;
					 fusible = false;
					 break;
			 }
			 is_fusible_[p] = fusible && fabs(xO-xFe)<1e-9;
			 // at present, cannot handle ferric or oxidized sulfur melts
		 }
	 }

	 // Compute end-member melt phase free energies
	 for (unsigned i=0; i<M_END; ++i)
	 {
		 Phase const &ph = ::phase[melt_endmember[i]];
		 Gfm_[i] = ph.model->Gf(ph, T_, P_);
	 }
 }

//-----------------------------------------------------------------------------//
/*! Calculate molar Gibbs free energy of melt.
 * 
 * \param x Mole fraction of each melt member. Energy returned will be for the
 * amount in x.
 */
double Melt_Model::Gfm(double const XM[E_END]) const
{
	double x[M_END];
	fill(x, x+M_END, 0.0);

	// Reorganize in analogy to CIPW norm into preferred melt phases.

	// 1: I have no phosphate melt

	// 2: I have no pyrite melt.

	// 3: I have no chromite melt.

	// 4: I hae no ilmenite melt.

	// 5: I have no fluorite melt.

	// 6: I have no calcite melt

	// 7: I have no zircon melt.

	// 8: Orthoclase

    double Q = 0.5*XM[E_K];
	x[M_KAlSi3O8] = Q;
	x[M_Al2O3] = 0.5*XM[E_AL] - Q;
	x[M_SiO2] = XM[E_SI] - 6*Q;

	// 9: Albite
	double Na2O = 0.5*XM[E_NA];
	Q = 0.5*min(Na2O, x[M_Al2O3]);
	x[M_NaAlSi3O8] = Q;
	Na2O -= 0.5*Q;
	x[M_Al2O3] -= 0.5*Q;
	x[M_SiO2] -= 3*Q;

	// 10: Anorthite
	x[M_CaO] = XM[E_CA];
	Q = min(x[M_CaO], x[M_Al2O3]);
	x[M_CaAl2Si2O8] = Q;
	x[M_CaO] -= Q;
	x[M_Al2O3] -= Q;
	x[M_SiO2] -= 2*Q;

	// 11: Acmite should come next, but I have no melt for it.

	// 12: Sodium metasilicate

	Q = Na2O;
	Na2O = 0.0;
	x[M_Na2SiO3] = Q;
	x[M_SiO2] -= Q;

	// 13: I have no melts for magnetite and am not presently implementing hematite.

	// 14-15: Magnesium diopside. I have no hedenbergite melt.

	x[M_MgO] = XM[E_MG];
	Q = min(x[M_MgO], x[M_CaO]);
	x[M_CaMgSi2O6] = Q;
	x[M_CaO] -= Q;
	x[M_MgO] -= Q;
	x[M_SiO2] -= 2*Q;

	// 16: Wollastonite

	Q = x[M_CaO];
	x[M_CaSiO3] = Q;
	x[M_CaO] = 0.0;
	x[M_SiO2] -= Q;

	// 17: Magnesium pyroxene (enstatite). I have no ferrosilite melt.

	Q = 0.5*x[M_MgO];
	x[M_Mg2Si2O6] = Q;
	x[M_MgO] = 0.0;
	x[M_SiO2] -= 2*Q;

	Q = 0.5*XM[E_FE];
	x[M_Fe2SiO4] = Q;
	x[M_SiO2] -= Q;

	if (x[M_SiO2]<0.0)
	{
		// Silica undersaturated

		double D = -x[M_SiO2];
		x[M_SiO2] = 0.0;

		// Enstatite to olivine.
		Q = min(D, x[M_Mg2Si2O6]);
		x[M_Mg2SiO4] = Q;
		x[M_Mg2Si2O6] -= Q;
		D -= Q;

		if (D>0.0)
		{
			// Albite to nepheline

			Q = min(0.5*D, x[M_NaAlSi3O8]);
			x[M_NaAlSiO4] += Q;
			x[M_NaAlSi3O8] -= Q;
			D -= 2*Q;

			if (D>0.0)
			{
				// Orthoclase to leucite

				Q = min(x[M_KAlSi3O8], D);
				x[M_KAlSi2O6] += Q;
				x[M_KAlSi3O8] -= Q;
				D -= Q;

				if (D>0.0)
				{
					// Do not have calcium orthosilicate. Go to lime instead.

					// Diopside to wollastonite and olivine.

					Q = min(x[M_CaMgSi2O6], 2*D);
					x[M_CaMgSi2O6] -= Q;
					x[M_CaSiO3] += Q;
					x[M_Mg2SiO4] += 0.5*Q;
					D -= 0.5*Q;


					if (D>0.0)
					{
						// Wollastonite to lime

						Q = min(D, x[M_CaSiO3]);
						x[M_CaO] += Q;
						x[M_CaSiO3] -= Q;
						D -= Q;

						if (D>0.0)
						{
							// Olivine to periclase

							Q = min(D, x[M_Mg2SiO4]);
							x[M_MgO] += 2*Q;
							x[M_Mg2SiO4] -= Q;
							D -= Q;

							// That's all I have melts for. If quartz is still deficient, this is not a valid melt.
							if (D<-1.0e-9)
							{
								return 1e10; // huge Gf turns off this melt.
							}
						}
					}
				}
			}
		}
	}

	x[M_H2O] = 0.5*XM[E_H];

//	cout << "Melt CIPW:" << endl;
//	for (unsigned m=0; m<M_END; ++m)
//	{
//		if (x[m]>0.0)
//		{
//			cout << phase[endmember[m]].name << ": " << x[m] << endl;
//		}
//	}

	// Compute melt free energy.

	double Gfm = 0.0;
	double Ntot = 0.0;
	for (unsigned i=0; i<M_END; ++i)
	{
		Ntot += mixN[i]*x[i];
	}
	for (unsigned i=0; i<M_END; ++i)
	{
		Gfm += Gfm_[i]*x[i];
		for (unsigned j=0; j<M_END; ++j)
		{
			Gfm += 0.5*x[i]*x[j]*W[i][j]/Ntot;
		}
	}

	// Entropy of mixing
	double nS = 0.0;
	for (unsigned i=0; i<M_END; ++i)
	{
		if (x[i]>0.0)
		{
			double const Nf = mixN[i]*x[i]/Ntot;
			nS += Ntot*Nlog(Nf);
        }
	}

	Gfm += R*T_*nS;

	return Gfm;
}

//-----------------------------------------------------------------------------//
double Melt_Model::Gf(double const XP[P_END]) const
{
	unsigned const NM = NP_;

	double xm[E_END];
	for (unsigned i=0; i<E_END; ++i)
	{
		xm[i] = Z_[i];
	}
	// Modify for amount of each fusible phase crystallized out
	double Result = Gfr_;
	for (unsigned i=0; i<NM; ++i)
	{
		double const x = XP[i];
		if (x>0.0)
		{
			Phase const &ph = phase_[ph_[i]];  
			unsigned const N = ph.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned z = ph.z[j];
				xm[z] -= x*ph.n[j];
				xm[z] = max(0.0, xm[z]);
			}
			Result += x*Gf_[ph_[i]];
		}
	}
    Result += this->Gfm(xm);

	return Result;
}

//-----------------------------------------------------------------------------//
/*! Compute the derivative of total free energy with the amount of a solid phase
 * extrated from the melt.
 * 
 * \param XP Number of moles of each solid phase. This is a compressed array of 
 * length Melt::NP_ specifying the number of moles of each active phase Melt::p_.
 * 
 * \param Gf0 Total free energy for the specified melt state, previously computed.
 * 
 * \param p Solid phase for which the derivative of total free energy is desired.
 * 
 * \return The derivate
 */
double Melt_Model::dGf(double XP[P_END], double const Gf0, unsigned const p) const
{
	unsigned const NM = NP_;

	// Compute current melt composition
	double xm[E_END];
	compute_current_melt_composition_(XP, xm);

	// See how much of the phase of interest can be crystallized from melt
	Phase const &ph = phase_[ph_[p]];  
	unsigned const N = ph.nz;
	double nmax = std::numeric_limits<double>::max();
	for (unsigned j=0; j<N; ++j)
	{
		unsigned z = ph.z[j];
		if (z != E_O)
	  	  nmax = min(nmax, xm[z]/ph.n[j]);
	}
	if (nmax>1.0e-9) // We can crystallize a significant quantity.
	{
		double old = XP[p];
		double h = min(nmax, 0.01);
		XP[p] += h;
		double Gf1 = Gf(XP);
		XP[p] = old;
		return (Gf1-Gf0)/h;
	}
	else if (XP[p]>1.0e-9) // We can melt a significant quantity.
	{
		double old = XP[p];
		double h = min(old, 0.01);
		XP[p] -= h;
		double Gf1 = Gf(XP);
		XP[p] = old;
		return (Gf0-Gf1)/h;
	}
	else // This phase wants to melt but is already all melted, or it wants
		// to crystallize but is already fully crystallized from melt. Change
		// of the phase is turned off so search direction reflects the constraint.
	{
		return 0.0;
	}
}

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
	unsigned const NPL = phase_.size();

	// Compute the total number of moles of atoms in the sample. This quantity will be
	// used in various tests of convergence.
	double cnorm = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		Check(Z_[i]>=0.0);
		cnorm += Z_[i];
	}
	
	double xm[E_END];
	compute_current_melt_composition_(XP, xm);

	DO:
	// Set temporary minimization set to all active phases.
	NP_ = NPL;
	for (unsigned i=0; i<NPL; ++i)
	{
		ph_[i] = i;
	}

	// Calculate free energy of initial guess of sample state, and its
	// gradient.

	double Gf = this->Gf(XP);
	cout << "Gf = " << fixed << setprecision(3) << Gf << endl;

	// Compute the gradient of the free energy with amount of each solid phase.
	// The gradient is a poor initial search direction, because it precipitates 
	// every solid phase that is supercooled, not just a set of phases that are 
	// mutually stable at the specified P and T. We therefore compute the
	// composition of the crystallizing material, and we will then perform
	// a solid phase minimization on this composition to obtain a set of such
	// mutually stable phases whose precipitation lowers the free energy of
	// the sample. 
	double xfreeze[E_END], p[NP_];
	fill(xfreeze, xfreeze+E_END, 0.0);
	cout << "Saturated extractable phases:" << endl;
	for (unsigned i=0; i<NP_; ++i)
	{
		double Gfp = is_fusible_[i]? dGf(XP, Gf, i) : 0.0;
		p[i] = -Gfp;
		Phase const &ph = phase_[i];

		// Accumulate conribution to freeze composition.
		if (p[i]>0.0)  // include only crystallizing phases
		{
 			unsigned n = ph.nz; 
			for (unsigned j=0; j<n; ++j)
			{
				unsigned z = ph.z[j];
				if (z != E_O && xm[z]<1.0e-9*cnorm)
				{
					cout << "  Phase "<< ph.name << " constrained from crystallizing" << endl;
					p[i] = 0.0;
					break;
				}
			}
			if (p[i]>0.0)
			{
				cout << ph.name << " " << p[i] << endl;
				for (unsigned j=0; j<n; ++j)
				{
					unsigned z = ph.z[j];
					xfreeze[z] += p[i]*ph.n[j];
				}
			}
		}
	}

	// Find the stable solid mineral assemblage for the crystallizing composition. 

	State state("1", T_, P_, xfreeze);
	auto status = state.do_ladder_update();

	// Report the phases in the crystallizing composition.
	cout << endl << "Constructing new search set. Crystallizing phases:" << endl;
	double const *const state_X = state.X();
	int const *const state_ph = state.ph();
	auto const &state_phase = state.phase();
	NP_ = 0;
	bool crystallizing[NPL];
	fill(crystallizing, crystallizing+NPL, false);
	double X[NPL];
	for (unsigned e=0; e<E_END; ++e)
	{
		int i = state_ph[e];
		if (state_X[e]>1.0e-9)
		{
			cout << " " << state_phase[state_ph[e]].name << " = " << state_X[e] << endl;
			ph_[NP_] = i;
			crystallizing[i] = true;
			X[NP_] = XP[i];
			NP_++;
		}
	}

	// Report the fusible phases already present as solids
	cout << "Melting phases: " << endl;
	for (unsigned i=0; i<NPL; ++i)
	{
		if (XP[i]>0.0 && !crystallizing[i])
		{
			ph_[NP_] = i;
			X[NP_] = XP[i];
			NP_++;
			cout << phase_[ph_[i]].name << endl;
		}
	}

	// Minimize on this melt set.

	double revised_Gf = minimize_trial_set_(X);
	compute_current_melt_composition_(X, xm);

	for (unsigned p=0; p<NPL; p++)
	{
		if (is_fusible_[p])
		{
			XP[p] = 0.0;
		}
	}
	for (unsigned p=0; p<NP_; ++p)
	{
		XP[ph_[p]] = X[p];
	}
	double mtot = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		if (xm[i]<1.0e-9*cnorm)
		{
			xm[i] = 0.0;
		}
		mtot += xm[i];
	}	
	if (mtot>0.0 && Gf - revised_Gf  > 1.0e-9*cnorm)
	{
		goto DO;
	}

	// Done. Construct new parent phase for melt.

	double V = 0.0;
	for (unsigned i=0; i<M_END; ++i)
	{
		double xi = xm[i];
		if (xi>0.0)
		{
			Phase const &ph = ::phase[melt_endmember[i]];
			V += ph.V*xi;
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

double Melt_Model::minimize_trial_set_(double p[P_END])
{
	unsigned const n = NP_;

	double g[n], h[n], xi[n];

	double xm[E_END];
	compute_current_melt_composition_(p, xm);

	double fp = Gf(p);
	double cnorm = accumulate(Z_, Z_+E_END, 0.0);
	double pnorm = 0.0;
	for (unsigned i=0; i<n; ++i)
	{
		xi[i] = dGf(p, fp, i);
		// Prune to enforce constraints
		if (xi[i]>0.0)
		{
			if (p[i]<=0.0)
			{
				// Prune melt of solid phase not present
				xi[i] = 0.0; 
			}
			else
			{
				cout << "Melting " << phase_[ph_[i]].name << endl;
			}
		}
		else if (xi[i]<0.0)
		{
			// Check for crystallization of element fully extracted
			Phase const &ph = phase_[ph_[i]];
			unsigned const Z = ph.nz;
			for (unsigned j=0; j<Z; ++j)
			{
				unsigned z = ph.z[j];
				if (z != E_O && xm[z]<1.0e-9*cnorm)
				{
					xi[i] = 0.0;
				}
			}
			if (xi[i] != 0.0)
			{
				cout << "Crystallizing " << ph.name << endl;
			}
		}
		pnorm += xi[i]*xi[i];
	}
	pnorm = sqrt(pnorm);

	for (int j=0; j<n; ++j)
	{
		g[j] = -xi[j];
		xi[j] = h[j] = g[j];
	}

	for (;;)
	{
		if (pnorm < 2.0e-7*cnorm) // To make sure search direction is meaningful
			break;

		// Compute freeze composition on pruned search direction. Clip upper
		// limit of search to not melt more of any phase than is actually present.
		// Also compute current melt composition.
		double xfreeze[E_END], xmelt[E_END];
		fill(xfreeze, xfreeze+E_END, 0.0);
		copy(Z_, Z_+E_END, xmelt);
		double x0 = 0.0;
		double x1 = numeric_limits<double>::max();
		for (unsigned i=0; i<n; ++i)
		{
			if (xi[i]!=0.0)
			{
				Phase const &ph = phase_[ph_[i]];
				unsigned Z = ph.nz;
				for (unsigned j=0; j<Z; ++j)
				{
					unsigned z = ph.z[j];
					xfreeze[z] += xi[i]*ph.n[j];
				}
			}
			if (p[i]!=0.0)
			{
				Phase const &ph = phase_[ph_[i]];
				unsigned Z = ph.nz;
				for (unsigned j=0; j<Z; ++j)
				{
					unsigned z = ph.z[j];
					xmelt[z] -= p[i]*ph.n[j];
				}
			}
			if (xi[i]<0.0)
			{
				x1 = min(x1, -p[i]/xi[i]);
			}
		}
		// Clip upper limit of search to not extract more freeze composition
		// than the amount of melt can supply.
		for (unsigned e=0; e<E_END; ++e)
		{
			if (e != E_O && xfreeze[e]>0.0)
			{
				x1 = min(x1, xmelt[e]/xfreeze[e]);
			}
		}

		double x2 = minimize(x0, x1, [&](double const e)
		                     {return Gfmelt(p,
		                                    xi,
		                                    e);});

		for (unsigned i=0; i<n; ++i)
		{
			p[i] += x2*xi[i];
			p[i] = max(0.0, p[i]);
		}

		fp = this->Gf(p);
	    compute_current_melt_composition_(p, xm);

		cout << endl;
		cout << "Line search completed. New Gf = " << fp << endl;
		cout << "Solid composition after line search:" << endl;
		for (unsigned i=0; i<n; ++i)
		{
			cout << "  " << phase_[ph_[i]].name << " = " << p[i] << endl;
		}
		cout << "Melt composition after line search:" << endl;
		for (unsigned i=0; i<E_END; ++i)
		{
			if (xm[i]>0.0)
			{
			cout << "  " << element_name[i] << " = " << xm[i] << endl;
			}
		}
		cout << endl;

		if (fabs(x2)*pnorm<3.0e-7*cnorm) 
			goto DONE;

		cout << "New pruned gradient calculation:" << endl;
		double dgg = 0., gg = 0.;
		for (unsigned i=0; i<n; ++i)
		{
			xi[i] = dGf(p, fp, i);

			// Prune to enforce constraints
			if (xi[i]>0.0)
			{
				if (p[i]<=0.0)
				{
					// Prune melt of solid phase not present
					xi[i] = 0.0; 
					cout << phase_[ph_[i]].name << " constrained from melting" << endl;
				}
				else
				{
					cout << "Melting " << phase_[ph_[i]].name << ", p = " << p[i] << endl;
				}
			}
			else if (xi[i]<0.0)
			{
				// Check for crystallization of element fully extracted
				Phase const &ph = phase_[ph_[i]];
				unsigned const Z = ph.nz;
				for (unsigned j=0; j<Z; ++j)
				{
					unsigned z = ph.z[j];
					if (z != E_O && xm[ph.z[j]]<1.0e-9*cnorm)
					{
						xi[i] = 0.0;
				     	cout << phase_[ph_[i]].name << " constrained from crystallizing" << endl;
					}
				}
				if (xi[i] != 0.0)
				{
					cout << "Crystallizing " << ph.name << ", p = " << p[i]  << endl;
				}
			}
			gg += g[i]*g[i];
			dgg += (xi[i]+g[i])*xi[i];
		}

		if (gg == 0.0)
			break;

		double  gam = dgg/gg;
		pnorm = 0.0;
		cout << endl;
		cout << "Conjugate gradient search direction:" << endl;
		for (unsigned i=0; i<n; i++)
		{
			g[i] = -xi[i];
			if (g[i] != 0.0)
			{
		   	   xi[i] = h[i] = g[i] + gam*h[i];
			}
			else
			{
		   	   xi[i] = h[i] = 0.0;
			}

			if (xi[i]!= 0.0)
			{
				cout << phase_[ph_[i]].name << " p = " << xi[i] << endl;
				pnorm += xi[i]*xi[i];
			}
		}
		pnorm = sqrt(pnorm);
	}
	
	DONE:
		return fp; 
#if 0
	{
		unsigned const NM = NP_;

		std::vector<double> x(M_END, 0.0);

		// Calculate composition of melt in terms of end members

		for (unsigned i=0; i<NM; ++i)
		{
			double const xi = min(Z_[i], max(0.0, p[i]));
			for (unsigned j=0; j<M_END; ++j)
			{
//				x[j] += xi*basis_[i][j];
			}
		}

		double sum = 0.0;
		for (unsigned i=0; i<M_END; ++i)
		{
			sum += x[i];
		}
		double rsum = 1.0/(sum + std::numeric_limits<double>::min());
		for (unsigned i=0; i<M_END; ++i)
		{
			x[i] *= rsum;
		}

		// Convert end member composition to elemental composition.

		double xe[E_END];
		fill(xe, xe+E_END, 0.0);
		double V = 0.0;
		for (unsigned i=0; i<M_END; ++i)
		{
			double xi = x[i];
			if (xi>0.0)
			{
				Phase const &ph = phase[melt_endmember[i]];
				unsigned nz = ph.nz;
				for (unsigned j=0; j<nz; ++j)
				{
					xe[ph.z[j]] += xi*ph.n[j];
				}
			    V += ph.V*xi;
			}
		}

		double Gfp = Gfm(x);
		
		Phase Result;
		Result.index = 0;
		Result.name = "mixed melt";
		Result.nz = 0;
		Result.V = V;
		for (unsigned i=0; i<E_END; ++i)
		{
			if (xe[i]>0.0)
			{
				Result.z[Result.nz] = i;
				Result.n[Result.nz] = xe[i];
				Result.nz++;
			}
		}

		Result.S0 = 0.0;
		Result.Hf0 = Gfp;
		Result.model = MELT;

		return Result;
	}
#endif
}

//-----------------------------------------------------------------------------//
double Melt_Model::Gfmelt(double const X[P_END], double p[M_END], double e) const
{
	unsigned const N = NP_;
	double x[N];
	for (unsigned i=0; i<N; ++i)
	{
		x[i] = X[i] + e*p[i];
	}
	return Gf(x);
}

//-----------------------------------------------------------------------------//
double State::melt(struct Phase &new_phase) const
{
	Melt_Model model(*this);

    // Initial guess is that all fusible phases are fully melted. 
	// We will then see what should crystallize out.

	double XP[P_END];
	fill(XP, XP+P_END, 0.0);
	
	new_phase = model.minimize_Gf(XP);

    double Gff = model.Gf(XP);

	return Gff; 
}


void Melt_Model::compute_current_melt_composition_(double const XP[], double xm[]) const
{
	for (unsigned i=0; i<E_END; ++i)
	{
		xm[i] = Z_[i];
	}
	for (unsigned i=0; i<NP_; ++i)
	{
		double const x = XP[i];
		if (x>0.0)
		{
			Phase const &ph = phase_[ph_[i]];  
			unsigned const N = ph.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned z = ph.z[j];
				xm[z] -= x*ph.n[j];
			}
		}
	}
}

char const * const endmember_element_name[] =
{
		"H",
		"Si",
		"Al",
		"Mg",
		"Fe(+2)",
		"Ca",
		"Na",
		"K"
};

