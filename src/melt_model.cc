// melt_model.cc
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
/*#include <algorithm>
#include <assert.h>
#include <fstream>
#include <vector>

#include "ds++/Assert.hh"
*/

#include "constants.hh"
#include "bm_melt.hh"
#include "phase.hh"

//inline double square(double x){ return x*x; }

using namespace std;

static const Melt myMelt;
extern Model const *const MELT = &myMelt;

//------------------------- Melts ---------------------------------------------

double Melt::Gf(Phase const &phase, double const T, double const P) const
{
	double const dT = T-T0;
	double const dT2 = T*T-T0*T0;
	double const drT = 1.0/T-1.0/T0;
	double const drT2 = 1.0/T/T - 1.0/T0/T0;
	double const sT = sqrt(T);
	double const sT0 = sqrt(T0);
	double const dsT = sT-sT0;
	double const drsT = 1.0/sT - 1.0/sT0;
	double const dlogT = log(T/T0);

	unsigned const N = phase.nz;

	double const Hf0 = phase.Hf0;
	double const S0 = phase.S0*1e-3; // to bring J to kJ
	double const V0 = phase.V;
	
	Melt_Phase const &sph = reinterpret_cast<Melt_Phase const &>(phase.data);

	double const A = sph.a;
	double const B = sph.b*1e-5;  // By convention, reported in units of 1e-5 kJ/K/K
	double const C = sph.c;
	double const D = sph.d;
	double const a0 = sph.a0*1e-5; // By convention, reported in units of 1e-5/K
	double const k0 = sph.k0;
	double const k0p = sph.k0p;
	double const k0pp = sph.k0pp;

	// Temperature terms
	double const Gt = Hf0 - T*S0 + A*dT + 0.5*B*dT2 - C*drT + 2*D*dsT
		- T*(A*dlogT + B*dT - 0.5*C*drT2 - 2*D*drsT);

	// Pressure term
	double const b = k0p*(2+k0p)/(k0*(1+k0p));
	double const c = 1.0/(k0p*(2+k0p));
	double const a = sqrt((1+c)/c);

	double const Vt = V0*(1 + a0*dT - 20*a0*dsT);

	double const Gfi =
		Gt + 
		P*Vt*(1-a+a*(1 - pow(1+b* P, 1-c))/(b*(c-1)* P))
		-P0*Vt*(1-a+a*(1 - pow(1+b*P0, 1-c))/(b*(c-1)*P0));


	return Gfi; 
}

double Melt::volume(Phase const &phase, double T, double P) const
{
  return phase.V;
}
