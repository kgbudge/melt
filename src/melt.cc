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
#include "melt.hh"
#include "Model.hh"
#include "phase.hh"
#include "State.hh"

//-----------------------------------------------------------------------------//
using namespace std;

//-----------------------------------------------------------------------------//
class Melt_Model
{
	public:
		Melt_Model(double const T, 
		           double const P, 
		           vector<Phase> const &phase,
		           vector<double> const &Gf,
		           State const &state) noexcept(false);

		unsigned nm() const noexcept { return nm_; }
		double x(unsigned i) const { Require(i<nm_); return x_[i]; }

		template<typename Real>
		Real Gfmelt(double XP[P_END]) const;

		template<typename Real>
		Real Gfm(Real const XM[E_END]) const;

		double Gfmelt(std::vector<double> const &X, double p[M_END], double e) const;
		
		Phase minimize_Gf(double XP[P_END]);

	private:
		template<typename Real>
		Real rebase_(double const XP[P_END], Real Gfm) const;
			 
		double T_; // temperature of sample. Pressure already folded into free energies.
		double x_[E_END]; // amount of each fusible element
		double Gf0_; // non-meltable phases total free energy
		unsigned nm_; // number of active phases
		unsigned p_[P_END]; // indices of active phases
		double Gfp_[P_END]; // molar free energies of active phases
		double Gfm_[M_END]; // molar free energies of endmember melts
};

//-----------------------------------------------------------------------------//
double const Wnorm = 1.0e-3;
double const Ww_q = 30967.*Wnorm;
double const Ww_cor = -16098.*Wnorm;
double const Ww_fo = 28874.*Wnorm;
double const Ww_fa = 35634.*Wnorm;
double const Ww_wo = 20375.*Wnorm;
double const Ww_sm = -96938.*Wnorm;
double const Ww_kal = 10374.*Wnorm;
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
double const Wq_kal = -33922.*Wnorm;
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
double const Wfo_fa = -37257.*Wnorm;
double const Wfo_wo = -31732.*Wnorm;
double const Wfo_sm = -41877.*Wnorm;
double const Wfo_kal = 22323.*Wnorm;
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
double const Wfa_wo = -12917.*Wnorm;
double const Wfa_sm = -90534.*Wnorm;
double const Wfa_kal = 23649.*Wnorm;
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
double const Wwo_sm = -13247.*Wnorm;
double const Wwo_kal = 17111.*Wnorm;
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
double const Wsm_kal = 6523.*Wnorm;
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
 * \param T Temperature (K)
 * \param P Pressure (kbar)
 * \param phase Current phase library
 * \param Gf Gibbs free energy of each phase at T and P (kJ/mol)
 * \param state Current state of the mineral ensemble.
 */ 
 Melt_Model::Melt_Model(double const T, 
                       double const P, 
                       vector<Phase> const &phase,
                       vector<double> const &Gf,
                       State const &state)
{
	// Start with the sample fully melted as our melt component. Determine
	// which fusible elements are present and compute a basis for expressing
	// the element composition of any fusible phase as primitive melt components.

	T_ = T;  // Temperature at which we're doing all this. Pressure is folded
	         // in to the free energies of phases and melt components already.
	 
	Gf0_ = 0.0;  // Gibbs free energy contribution of nonfusible phases 
    fill(x_, x_+E_END, 0.0); // To accumulate melt composution
	 
	cout << "Sample phases:" << endl;
	for (unsigned i=0; i<E_END; ++i)
	{
		double const x = state.x[i];
		double xpr[E_END];  // to store composition of phase
		fill(xpr, xpr+E_END, 0.0); 
		if (x>0.0)  // Is this phase actually present?
		{
			unsigned p = state.p[i];
			Phase const &ph = phase[p];	
			cout << ph.name << endl;
			unsigned const N = ph.nz; // Number of elements in the phase
			double xO = 0.0;          // To accumulate oxygen balance of the phase
			bool fusible = true;
			for (unsigned j=0; j<N && fusible; ++j) 
			{
				double const xj = ph.n[j];  // Number of moles of element j in the phase.
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
				Gf0_ += state.x[i]*Gf[p];
			}
			else
			{
				for (unsigned m=0; m<E_END; ++m)
				{
					x_[m] += xpr[m];
				}
			}
		}
	}
	cout << "Full melt elemental molar composition:" << endl << defaultfloat;
	for (unsigned i=0; i<E_END; ++i)
	{
		if (x_[i]>0)
	  	  cout << element_name[i] << ": " << x_[i] << endl;
	}

	 // Now go over all phases that might interact with melt, and add them to the
	 // phase list.
	cout << "Active phases:" << endl;
    nm_ = 0; // accumulate number of active phases
	for (unsigned i=0; i<P_END; ++i)
	{
		Phase const &ph = phase[i];	
		unsigned const N = ph.nz; // Number of elements in the phase
		double xO = 0.0, xFe = 0.0;          // To accumulate oxygen balance of the phase
		bool fusible = (N>0);
		for (unsigned j=0; j<N && fusible; ++j) 
		{
			double const xj = ph.n[j];  // Number of moles of element j in the phase.
			if (ph.z[j] != E_O && xj>0.0 && x_[ph.z[j]] == 0.)
			{
				fusible = false;
				break;
			}
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
		}
		if (fusible && fabs(xO-xFe)<=1e-9) // at present, cannot handle ferric or oxidized sulfur melts
		{
			cout << phase[i].name << endl;
			p_[nm_] = i;
			Gfp_[nm_] = Gf[i];
			nm_++;
		}
	}

	 // Melt phase molar free energies
	 for (unsigned i=0; i<M_END; ++i)
	 {
		 Gfm_[i] = Gf[endmember[i]];
	 }
}

//-----------------------------------------------------------------------------//
/*! Calculate molar Gibbs free energy of melt.
 * 
 * \param x Mole fraction of each melt member. Energy returned will be for the
 * amount in x.
 */
template<typename Real>
Real Melt_Model::Gfm(Real const XM[E_END]) const
{
	unsigned const NM = nm_;

	Real zero = to_Real<Real>(0.0, E_END);

	Real x[M_END];
	fill(x, x+M_END, zero);

	// Reorganize in analogy to CIPW norm into preferred melt phases.

	// 1: I have no phosphate melt

	// 2: I have no pyrite melt.

	// 3: I have no chromite melt.

	// 4: I hae no ilmenite melt.

	// 5: I have no fluorite melt.

	// 6: I have no calcite melt

	// 7: I have no zircon melt.

	// 8: Orthoclase

    Real Q = 0.5*XM[E_K];
	x[M_KAlSi3O8] = Q;
	x[M_Al2O3] = 0.5*XM[E_AL] - Q;
	x[M_SiO2] = XM[E_SI] - 6*Q;

	// 9: Albite
	Real Na2O = 0.5*XM[E_NA];
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
	Na2O = zero;
	x[M_Na2SiO3] = Q;
	x[M_SiO2] -= Q;

	// 13: I have no melts for magnetite and am not presently implementing hematite.

	// 14-15: Magnesium diopside. I have no hedenbergite melt.

	x[M_MgO] = XM[E_MG];
	Q = min(x[M_MgO], XM[M_CaO]);
	x[M_CaMgSi2O6] = Q;
	x[M_CaO] -= Q;
	x[M_MgO] -= Q;
	x[M_SiO2] -= 2*Q;

	// 16: Wollastonite

	Q = x[M_CaO];
	x[M_CaSiO3] = Q;
	x[M_CaO] = zero;
	x[M_SiO2] -= Q;

	// 17: Magnesium pyroxene (enstatite). I have no ferrosilite melt.

	Q = 0.5*x[M_MgO];
	x[M_Mg2Si2O6] = Q;
	x[M_MgO] = zero;
	x[M_SiO2] -= 2*Q;

	Q = 0.5*XM[E_FE];
	x[M_Fe2SiO4] = Q;
	x[M_SiO2] -= Q;

	if (x[M_SiO2]<0.0)
	{
		// Silica undersaturated

		Real D = -x[M_SiO2];
		x[M_SiO2] = zero;

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
								return to_Real<Real>(1e10, NM); // huge Gf turns off this melt.
							}
						}
					}
				}
			}
		}
	}

	x[M_H2O] = 0.5*XM[E_H];

	cout << "Melt CIPW:" << endl;
	for (unsigned m=0; m<M_END; ++m)
	{
		if (x[m]>0.0)
		{
			cout << phase[endmember[m]].name << ": " << (double)x[m] << endl;
			for (unsigned i=0; i<E_END; ++i)
			{
				cout << "dx[" << i << "] = " << dydx(x[m], i) << endl;
			}
		}
	}

	// Compute melt free energy.

	Real Gfm = zero;
	Real Ntot = zero;
	for (unsigned i=0; i<M_END; ++i)
	{
		Ntot += mixN[i]*x[i];
	}
	double const T = T_;
	for (unsigned i=0; i<M_END; ++i)
	{
		Gfm += Gfm_[i]*x[i];
		for (unsigned j=0; j<M_END; ++j)
		{
			Gfm += 0.5*x[i]*x[j]*W[i][j]/Ntot;
		}
	}
	
	cout << "Gfm (before entropy): " << value(Gfm) << endl;
	for (unsigned i=0; i<E_END; ++i)
	{
		cout << "dGfm[" << i << "] = " << dydx(Gfm, i) << endl;
	}

	// Entropy of mixing
	Real nS = zero;
	for (unsigned i=0; i<M_END; ++i)
	{
		if (x[i]>0.0)
		{
			Real const Nf = mixN[i]*x[i]/Ntot;
			nS += Ntot*Nlog(Nf);
        }
	}

	Gfm += R*T*nS;
	
	cout << "Gfm (after entropy): " << value(Gfm) << endl;
	for (unsigned i=0; i<E_END; ++i)
	{
		cout << "dGfm[" << i << "] = " << dydx(Gfm, i) << endl;
	}


	return Gfm;
}

//-----------------------------------------------------------------------------//
template<typename Real>
Real Melt_Model::Gfmelt(double XP[P_END]) const
{
	unsigned const NM = nm_;

	Real xm[E_END];
	for (unsigned i=0; i<E_END; ++i)
	{
		xm[i] = to_Real<Real>(x_[i], i, E_END);
	}
	// Modify for amount of each fusible phase crystallized out
	for (unsigned i=0; i<NM; ++i)
	{
		if (XP[i]>0.0)
		{
			Phase const &ph = phase[p_[i]];  
			unsigned const N = ph.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned z = ph.z[j];
				xm[z] -= ph.n[j]*XP[i];
			}
		}
	}
    Real Gfm = this->Gfm(xm);

	// Gfm now contains the total free energy of the melt plus its derivatives with
	// molar composition by cation.

	// Next we compute the total free energy of the sample (melt plus phases) and
	// the derivative of the free energy of the sample with each phase
	// crystallizing out of it.

	return rebase_<Real>(XP, Gfm);
}

template<>
double Melt_Model::rebase_(double const XP[P_END], double Gfm) const
{
	unsigned const NM = nm_;

	double Gft = Gf0_ + Gfm;  // Inactive phases contribution to free energy

	for (unsigned i=0; i<NM; ++i)
	{
		double const xi = min(x_[i], max(0.0, XP[i]));
		Gft += xi*Gfp_[i];
	}

	return Gft;
}

template<>
D1 Melt_Model::rebase_(double const XP[P_END], D1 Gfm) const
{
	unsigned const NM = nm_;

	D1 zero(0.0, NM);
	D1 Gft(Gf0_, NM);  // Inactive phases contribution to free energy

	Gft += value(Gfm); // will add derivatives later
	for (unsigned i=0; i<NM; ++i)
	{
		D1 const xi(min(x_[i], max(0.0, XP[i])), i, NM);
		Gft += xi*Gfp_[i];

		// Now for derivatives
		Phase const &ph = phase[p_[i]];
		unsigned const N = ph.nz;
		double deriv = 0;
		for (unsigned j=0; j<N; ++j)
		{
			unsigned z = ph.z[j];
		    deriv -= ph.n[j]*Gfm.dydx(z);
		}
		Gft.dydx(i, deriv + Gfp_[i]);
	}

	cout << "Gft: " << value(Gft) << endl;
	for (unsigned i=0; i<NM; ++i)
	{
		cout << "dGfs[" << i << "] = " << dydx(Gft, i) << endl;
	}

	// Now do it the hard way
	double dGft = value(Gft);
	double Xp[P_END];
	copy(XP, XP+NM, Xp);
	for (unsigned i=0; i<NM; ++i)
	{
		double old = Xp[i];
		Xp[i] += 0.01;
		double f2 = Gfmelt<double>(Xp);
		cout << "dGfs[" << i << "] = " << dydx(Gft, i) << endl;
		cout << "    [" <<  i << "] = " << (f2-dGft)*100. << endl;
		Xp[i] = old;
	}

	// Now compute total free energy. 
	return Gft;
}

//-----------------------------------------------------------------------------//
/*! Find the degree of melting of each solid phase that minimized total Gibbs free energy.
 * 
 * \param[in,out] X Contains initial guess of amount of melt of each phase. On return,
 * contains amount of melt of each phase that minimizes free energy.
 */
Phase Melt_Model::minimize_Gf(double XP[P_END])
{
	  // Compute the initial gradient and Hessian matrix

	unsigned const N = nm_;

	double g[N], xi[N], h[N], p[N];

		// Calculate initial Gf and gradient.

	D1 Gf2 = Gfmelt<D1>(XP);

#if 0
    double Gf = Gf2.y();

	cout << "Melt parameters:" << endl;
	cout << "  Gf = " << Gf << endl;
	for (unsigned i=0; i<N; ++i)
	{
		cout << "  X[" << i << "] = " << X[i] << endl;
		p[i] = xi[i] = h[i] = g[i] = -Gf2.dydx(i);
		cout << "dDf[" << i << "] = " << -xi[i] << endl;
	}

	DO:
	double x0 = -numeric_limits<double>::max();
	double x1 = numeric_limits<double>::max();
	double norm = 0.0;
	double cnorm = 0.0;
	for (unsigned i=0; i<N; ++i)
	{
		// Prune value up against a compositional constraint
		if (p[i]>0.0 && fabs(X[i]-x_[i])<1.0e-9)
		{
			X[i] = x_[i];
			p[i] = 0.0;
		}
		if (p[i]<0.0 && fabs(X[i])<1.0e-9)
		{
			X[i] = 0.0;
			p[i] = 0.0;
		}
		if (p[i]!=0.0)
		{
			double e = -X[i]/p[i];
			if (p[i]<0) x1 = min(x1, e); else x0 = max(x0, e);
			e = (x_[i]-X[i])/p[i];
			if (p[i]>0) x1 = min(x1, e); else x0 = max(x0, e);
		}
		norm += p[i]*p[i];
		cnorm += x_[i]*x_[i];
	}
	norm = sqrt(norm);
	cnorm = sqrt(cnorm);

	if (norm > 2.0e-7*cnorm) // To make sure gradient is meaningful
	{
		double x2 = minimize(x0, x1, [&](double const e)
			                     {return Gfmelt(X,
			                                    p,
			                                    e);});

		if (fabs(x2)*norm>3.0e-7*cnorm)
		{
			double sum = 0.0;
			for (unsigned i=0; i<N; ++i)
			{
				X[i] += x2*p[i];
				X[i] = min(x_[i], max(0.0, X[i]));
			}

			Gf2 = Gfmelt<D1>(X);

			cout << "Melt parameters:" << endl;
			cout << "  Gf = " << Gf << endl;
			double dgg = 0., gg = 0.;
			for (unsigned i=0; i<N; ++i)
			{
				cout << "  X[" << i << "] = " << X[i] << endl;
				xi[i] = Gf2.dydx(i);
				cout << "dDf[" << i << "] = " << xi[i] << endl;
				gg += g[i]*g[i];
				dgg += (xi[i]+g[i])*xi[i];
			}

			if (gg == 0.0)
				goto DONE;

			double  gam = dgg/gg;
			for (unsigned i=0; i<N; i++)
			{
				g[i] = -xi[i];
				p[i] = xi[i] = h[i] = g[i] + gam*h[i];
			}
			
			goto DO;
		}
	}

	DONE:
	{
		unsigned const NM = nm_;

		std::vector<double> x(M_END, 0.0);

		// Calculate composition of melt in terms of end members

		for (unsigned i=0; i<NM; ++i)
		{
			double const xi = min(x_[i], max(0.0, X[i]));
			for (unsigned j=0; j<M_END; ++j)
			{
				x[j] += xi*basis_[i][j];
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
				Phase const &ph = phase[endmember[i]];
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
double Melt_Model::Gfmelt(std::vector<double> const &X, double p[M_END], double e) const
{
#if 0
	unsigned const N = nm_;
	std::vector<double> x(N);
	for (unsigned i=0; i<N; ++i)
	{
		x[i] = X[i] + e*p[i];
	}
	return Gfmelt<double>(x);
#endif
}

//-----------------------------------------------------------------------------//
double melt(double const T, 
            double const P, 
            vector<Phase> const &phase,
            vector<double> const &Gf,
            State const &state, 
            struct Phase &new_phase)
{
	Melt_Model model(T, P, phase, Gf, state);

    // Initial guess is that all fusible phases are fully melted. 
	// We will then see what should crystallize out.

	double XP[P_END];
	fill(XP, XP+P_END, 0.0);
	
	new_phase = model.minimize_Gf(XP);

	double Gff = model.Gfmelt<double>(XP);

	return Gff; 
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

