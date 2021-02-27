/*
 * Melt_Model__Gfm.cc
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

#include <cmath>

#include "constants.hh"

//-----------------------------------------------------------------------------//
using namespace std;

//-----------------------------------------------------------------------------//
double const Wnorm = 1.0e-3;
double const Ww_q = 0.*Wnorm;  // Melting point 1713C. Ours: 1719C
double const Ww_cor = -0.*Wnorm;  // Melting point 2044C.  Ours: 2042C
double const Ww_fo = 0.*Wnorm;  // Melting point 1890C. Ours: 1823C
double const Ww_fa = 0.*Wnorm; // Incongruently melts to liquid and solid iron at 1217 C. 
                               // Ours: 1195C to liquid and wustite.
double const Ww_wo = 0.*Wnorm; // Wikipedia: melts at 1540C Ours: 1558C
double const Ww_sm = -0.*Wnorm;  // Wikipedia: 1088C Ours: 1087C
double const Ww_kal = 0.*Wnorm;  // 1686C Ours: 1848C
double const Ww_co2 = 0*Wnorm;
double const Ww_ne = 0*Wnorm; // 1526C Ours: 1527C
double const Ww_ab = 0*Wnorm; // 1115C Ours: 1116C incongruent
double const Ww_pc = 0*Wnorm;
double const Ww_s = 0*Wnorm;
double const Ww_ha = 0*Wnorm;
double const Ww_en = 0*Wnorm;
double const Ww_li = 0*Wnorm;
double const Ww_mdi = 0*Wnorm;
double const Ww_an = 0*Wnorm;
double const Ww_or = 0*Wnorm;
double const Wq_cor = -2520*Wnorm; // calibrated from SiO2-Al2O3 eutectic at 1868K=1595C  Al 0.9m vs. 0.4m
double const Wq_fo = 0; // do not coexist
double const Wq_fa = 2880.*Wnorm; // calibrated from Q-Fa eutectic near 1422K=1149C
double const Wq_wo = 1850.*Wnorm;// calibrated from eutectic of 1699K=1426C
double const Wq_sm = -3050*Wnorm;// calibrated from eutectic at 1062K=789C 
double const Wq_kal = -0.*Wnorm; // do not coexist
double const Wq_co2 = 0*Wnorm; 
double const Wq_ne = 0*Wnorm; // do not coexist

double const Wq_ab = -6660*Wnorm;  // eutectic at 1062 C
double const Wsm_ab = -1000*Wnorm;

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

double Nlog(double const &x)
{
	if (x>0.0)
	{
		return x*log(x);
	}
	else
	{
		return 0.0;
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
	compute_melt_endmember_composition_(XM, x);

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
