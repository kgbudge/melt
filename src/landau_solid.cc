// landau_solid.cc
//
// Copyright (C) 2021 Kent G. Budge
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

#include "landau_solid.hh"

#include <cmath>
#include <algorithm>

#include "constants.hh"
#include "Model.hh"
#include "Phase.hh"

inline double square(double x){ return x*x; }

using namespace std;

static Landau_Solid mySolid;
extern Model const *const LANDAU_SOLID = &mySolid;

//------------------------- Solids ---------------------------------------------

double Landau_Solid::Gf(Phase const &phase, double const T, double const P) const
{
	// Note that Gf is the apparent free energy of formation (in kJ/mol),
	// which ignores entropy of the elements.

	// Temperature is in K and P in kbar.

	double Gfi = Solid::Gf(phase, T, P);

	// Disorder is a perturbation on a regular solid.

	Solid_Phase const &sph = reinterpret_cast<Solid_Phase const &>(phase.data);

	double Tc = sph.Tc;
	double Smax = sph.Smax*1e-3;  // J to kJ
	double Vmax = sph.Vmax;

	double Tcs = Tc + Vmax*P/Smax;
	double Tcs0 = Tc + Vmax*P0/Smax;

	double Q0 = sqrt(sqrt((Tcs0-T0)/Tc));
	double Q = sqrt(sqrt(max(0.0, Tcs-T)/Tc));

	double Gdis = Tc*Smax*(Q0*Q0 - Q0*Q0*Q0*Q0*Q0*Q0/3.) - Smax*(Tcs*Q*Q-Tc*Q*Q*Q*Q*Q*Q/3.)
		- T*(Smax*(Q0*Q0-Q*Q)) + P*(Vmax*Q0*Q0);

	Gfi += Gdis;
	
	return Gfi;
}

double Landau_Solid::volume(Phase const &phase, double T, double P) const
{
  return phase.V;
}
