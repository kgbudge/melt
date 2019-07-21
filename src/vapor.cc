// vapor.cc
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
//#include <algorithm>
//#include <assert.h>
//#include <fstream>
//#include <vector>

//#include "ds++/Assert.hh"

#include "constants.hh"
#include "phase.hh"
#include "vapor.hh"

inline double square(double x){ return x*x; }

static Vapor const myVapor;
Model const *const VAPOR = &myVapor;

using namespace std;

//------------------------- Vapor ----------------------------------------------

double Vapor::Gf(Phase const &phase, double const T, double const P) const
{
	// Note that Gf is the apparent free energy of formation (in kJ/mol),
	// which ignores entropy of the elements.

	// Temperature is in K and P in kbar.

	double const S0 = phase.S0 * 1e-3;
	double const Hf0 = phase.Hf0;

	Vapor_Phase const &sph = reinterpret_cast<Vapor_Phase const &>(phase.data);

	double const a = sph.a;
	double const b = sph.b * 1.0e-5;
	double const c = sph.c;
	double const d = sph.d;

	double const Gt = Hf0 - T*S0 + a*(T-T0) + 0.5*b*(T*T-T0*T0) - c*(1/T-1/T0) + 2*d*(sqrt(T)-sqrt(T0))
		- T*(a*log(T/T0) - b*(T-T0) + 0.5*c*(1/(T*T)-1/(T0*T0)) + 2*d*(1/sqrt(T)-1/sqrt(T0))); 

	double const Gf = Gt + R*T*log(P/P0);
	return Gf;
}

double Vapor::volume(Phase const &phase, double const T, double const P) const
{
	return  R*T/Vapor::P(phase, T, P);
}

double Vapor::P(Phase const &phase, double const T, double const Gf) const
{
	double const S0 = phase.S0 * 1e-3;
	double const Hf0 = phase.Hf0;

	Vapor_Phase const &sph = reinterpret_cast<Vapor_Phase const &>(phase.data);
	
	double const a = sph.a;
	double const b = sph.b * 1.0e-5;
	double const c = sph.c;
	double const d = sph.d;
	
	double const Gt = Hf0 - T*S0 + a*(T-T0) + 0.5*b*(T*T-T0*T0) - c*(1/T-1/T0) + 2*d*(sqrt(T)-sqrt(T0))
		- a*log(T/T0) + b*(T-T0) - 0.5*c*(1/(T*T)-1/(T0*T0)) - 2*d*(1/sqrt(T)-1/sqrt(T0)); 

	double const P = P0*exp((Gf - Gt)/(R*T));
	return P;
}
