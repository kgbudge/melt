// solid.cc
//
// Copyright (C) 20169- Kent G. Budge
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

#include "ds++/Assert.hh"*/

#include "constants.hh"
#include "model.hh"
#include "phase.hh"
#include "solid.hh"

inline double square(double x){ return x*x; }

using namespace std;

static Solid mySolid;
extern Model const *const SOLID = &mySolid;

//------------------------- Solids ---------------------------------------------

double Solid::Gf(Phase const &phase, double const T, double const P) const
{
	// Note that Gf is the apparent free energy of formation (in kJ/mol),
	// which ignores entropy of the elements.

	// Temperature is in K and P in kbar.

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

	Solid_Phase const &sph = reinterpret_cast<Solid_Phase const &>(phase.data);
	
	double const A = sph.a;
	double const B = sph.b*1e-5;  // By convention, reported in units of 1e-5 kJ/K/K
	double const C = sph.c;
	double const D = sph.d;
	double const a0 = sph.a0*1e-5; // By convention, reported in units of 1e-5/K
	double const k0 = sph.k0;
	double const k0p = sph.k0p;
	double const k0pp = sph.k0pp;

	if (sph.Tlow>0.0 && 
	    (T<sph.Tlow || T>sph.Thigh || P<sph.Plow || P>sph.Phigh))
	{
		return 1e5;
	}
	
	// Temperature terms
	double const Gt = Hf0 + A*dT + 0.5*B*dT2 - C*drT + 2*D*dsT
		          - T*(S0 + A*dlogT + B*dT - 0.5*C*drT2 - 2*D*drsT);

	// Pressure term
	double const b = k0p*(2+k0p)/(k0*(1+k0p));
	double const c = 1.0/(k0p*(2+k0p));
	double const a = sqrt((1+c)/c);

	unsigned n = 0;
	for (unsigned j=0; j<N; ++j)
	{
		n += phase.n[j];
	}

	double const theta = 10636/(1e3*S0/n+6.44);

	double const xi0 = square(theta/(T0*expm1(theta/T0)))*exp(theta/T0);

	double const Pth = a0*k0*theta*(1/expm1(theta/T)-1/expm1(theta/T0))/xi0;

	double Vi = 1-a*(1-pow(1+b*(P0-Pth), -c));
	double Vf = 1-a*(1-pow(1+b*(P-Pth), -c));
	if (1 >= b*Pth && 1 >= -b*(P-Pth) && Vi>0 && Vf>0)
	{
		double const Gfi = 
			Gt + 
			P*V0*(1-a+a*(pow(1-b*Pth, 1-c) - pow(1+b*(P-Pth), 1-c))/(b*(c-1)*P))
			-P0*V0*(1-a+a*(pow(1-b*(Pth), 1-c) - pow(1+b*(P0-Pth), 1-c))/(b*(c-1)*P0));

		return Gfi;
	}
	else
	{
		return 1.0e5; // off range of validity of table
	}
}

double Solid::volume(Phase const &phase, double T, double P) const
{
  return phase.V;
}
