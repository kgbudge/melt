// bm_solid.cc
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
#include "bm_solid.hh"
#include "constants.hh"
#include "Phase.hh"

static BM_Solid const myBM_Solid;
Model const *const BM_SOLID = &myBM_Solid;

//------------------------- Birch-Murnaghan Melts ---------------------------------------------

double BM_Solid::Gf(Phase const &phase, double const T, double const P) const
{
	unsigned const N = phase.nz;

	double const dH_sTrPr = phase.Hf0;
    double const S_sTrPr = phase.S0*1e-3; // J to kJ 
    double const Vr = phase.V;
	
	BM_Solid_Phase const &sph = reinterpret_cast<BM_Solid_Phase const &>(phase.data);

	double c0 = sph.c0;
	double c1 = sph.c1;
	double c2 = sph.c2;
	double c3 = sph.c3;
	double dvdt = sph.dvdt;
	double dvdp = sph.dvdp;
	double d2vdpdt = sph.d2vdpdt;
	double Kp = sph.Kp;

	double const Tr = T0;
	double const Pr = P0;

	// Integrate the molar enthalpy and entropy of the solid from Tr, Pr to T, Pr

	double const dH_TPr = dH_sTrPr + c0*(T-Tr) + 2*c1*(sqrt(T)-sqrt(Tr))
		- c2*(1/T - 1/Tr) - 0.5*c3*(1/(T*T)-1/(Tr*Tr));
	
	double const S_TPr = S_sTrPr + c0*(log(T)-log(Tr)) - 2*c1*(1/sqrt(T)-1./sqrt(Tr)) 
		- 0.5*c2*(1/(T*T) - 1/(Tr*Tr)) - (1/3.)*c3*(1/(T*T*T)-1/(Tr*Tr*Tr));

	// giving the Gibbs free energy:

	double const dG_TPr = dH_TPr - T*S_TPr; 

	// Now the pressure integration

	double const V0 = Vr + dvdt*(T-Tr);
	double Vl = V0;

	p_phase = this;
	p_T = T;
	p_V = Vr;

	solve(P, Vl, [=](double V){
		return p(sph, V);});

	double const K = -V0/(dvdp+d2vdpdt*(T-Tr));

	double const rdl = pow(V0/Vl, 2./3.) - 1.;
	double const dG_TP  = dG_TPr  + P*Vl - Pr*V0 + (9./8.)*K*V0*rdl*rdl*(1+0.5*(Kp-4)*rdl);

	return dG_TP; 
}

double BM_Solid::volume(Phase const &phase, double T, double P) const
{
  return phase.V;
}

double BM_Solid::p_T, BM_Solid::p_V;
BM_Solid const *BM_Solid::p_phase;

double BM_Solid::p(BM_Solid_Phase const &sph, double Vl)
{
	double c0 = sph.c0;
	double c1 = sph.c1;
	double c2 = sph.c2;
	double c3 = sph.c3;
	double dvdt = sph.dvdt;
	double dvdp = sph.dvdp;
	double d2vdpdt = sph.d2vdpdt;
	double Kp = sph.Kp;

	double const T = p_T;
	double const Tr = T0;
	double const Vr = p_phase->p_V;

	double const V0 = Vr + dvdt*(T-Tr);
  double const K = -V0/(dvdp+d2vdpdt*(T-Tr));
  return 1.5*K*(pow(V0/Vl, 7./3.)-pow(V0/Vl, 5./3.))*(1 - 0.75*(4-Kp)*(pow(V0/Vl, 2./3.)-1)); 
}

