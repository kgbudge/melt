// bm_model.cc
//
// Copyright (C) 2019 - Kent G. Budge
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <cmath>

#include "algorithm.hh"
#include "bm_melt.hh"
#include "constants.hh"
#include "phase.hh"

static BM_Melt const myBM_Melt;
Model const *const BM_MELT = &myBM_Melt;

//------------------------- Birch-Murnaghan Melts ---------------------------------------------

double BM_Melt::Gf(Phase const &phase, double const T, double const P) const
{
	unsigned const N = phase.nz;

	double const dH_sTrPr = phase.Hf0;
    double const S_sTrPr = phase.S0*1e-3; // J to kJ 
    double const Vr = phase.V;
	
	Melt_Phase const &sph = reinterpret_cast<Melt_Phase const &>(phase.data);

	double c0 = sph.c0;
	double c1 = sph.c1;
	double c2 = sph.c2;
	double c3 = sph.c3;
	double Tf = sph.Tf;
	double dS_f = sph.dS_f;
	double Cpl = sph.Cpl;
	double dvdt = sph.dvdt;
	double dvdp = sph.dvdp;
	double d2vdpdt = sph.d2vdpdt;
	double Kp = sph.Kp;

	double const Tr = T0;
	double const Pr = P0;

	// Integrate the molar enthalpy and entropy of the solid from Tr, Pr to Tf, Pr

	double const dH_sTfPr = dH_sTrPr + c0*(Tf-Tr) + 2*c1*(sqrt(Tf)-sqrt(Tr)) - c2*(1/Tf - 1/Tr) - 0.5*c3*(1/(Tf*Tf)-1/(Tr*Tr));
	double const S_sTfPr = S_sTrPr + c0*(log(Tf)-log(Tr)) - 2*c1*(1/sqrt(Tf)-1./sqrt(Tr)) - 0.5*c2*(1/(Tf*Tf) - 1/(Tr*Tr)) - (1/3.)*c3*(1/(Tf*Tf*Tf)-1/(Tr*Tr*Tr));

	// Now from solid to melt, with constraint that dG is zero (because melting is reversible)

	double const dH_lTfPr = dH_sTfPr + Tf*dS_f;
	double const S_lTfPr = S_sTfPr + dS_f;

	// Now from melt at Tf, Pr to melt at T, Pr

	double const dH_lTPr = dH_lTfPr + Cpl*(T-Tf);
	double const S_lTPr = S_lTfPr + Cpl*(log(T)-log(Tf));

	// giving the Gibbs free energy:

	double const dG_lTPr = dH_lTPr - T*S_lTPr; 

	// Now the pressure integration

	double const V0 = Vr + dvdt*(T-Tr);
	double Vl = V0;

	p_phase = this;
	p_T = T;
	p_V = Vr;

	solve(P, Vl, [=](double V){
		return p(sph, V);});

	double const K = V0/(dvdp+d2vdpdt*(T-Tr));

	double const rdl = pow(V0/Vl, 2./3.) - 1.;
	double const dG_lTP = dG_lTPr + P*Vl - Pr*V0 + (9./8.)*K*V0*rdl*rdl*(1+0.5*(Kp-4)*rdl);

	return dG_lTP; 
}

double BM_Melt::volume(Phase const &phase, double T, double P) const
{
  return phase.V;
}

double BM_Melt::p_T, BM_Melt::p_V;
BM_Melt const *BM_Melt::p_phase;

double BM_Melt::p(Melt_Phase const &sph, double Vl)
{
	double c0 = sph.c0;
	double c1 = sph.c1;
	double c2 = sph.c2;
	double c3 = sph.c3;
	double Tf = sph.Tf;
	double dS_f = sph.dS_f;
	double Cpl = sph.Cpl;
	double dvdt = sph.dvdt;
	double dvdp = sph.dvdp;
	double d2vdpdt = sph.d2vdpdt;
	double Kp = sph.Kp;

	double const T = p_T;
	double const Tr = T0;
	double const Vr = p_phase->p_V;

	double const V0 = Vr + dvdt*(T-Tr);
  double const K = V0/(dvdp+d2vdpdt*(T-Tr));
  return 1.5*K*(pow(V0/Vl, 7./3.)-pow(V0/Vl, 5./3.))*(1 - 0.75*(4-Kp)*(pow(V0/Vl, 2./3.)-1)); 
}

