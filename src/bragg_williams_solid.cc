// bragg_williams_solid.cc
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

#include "algorithm.hh"
#include "constants.hh"
#include "model.hh"
#include "phase.hh"
#include "bragg_williams_solid.hh"

inline double square(double x){ return x*x; }

using namespace std;

static Bragg_Williams_Solid mySolid;
extern Model const *const BRAGG_WILLIAMS_SOLID = &mySolid;

//------------------------- Solids ---------------------------------------------

double Bragg_Williams_Solid::Gf(Phase const &phase, double const T, double const P) const
{
	// Temperature is in K and P in kbar.

	double Gfi = Solid::Gf(phase, T, P);
	// Disorder is a perturbation on a regular solid.

	Solid_Phase const &sph = reinterpret_cast<Solid_Phase const &>(phase.data);

	double dH = sph.dH;
	double dV = sph.dV;
	double Wh = sph.W;
	double Wv = sph.Wv;
	double n = sph.n;
	double Fac = sph.Fac;

	double W = Wh+P*Wv;

	double Tc = 2*W/(R*(1+n));

	double a = 0.0, b = 1.0;
	
	solvebr(a, b, [=](double Q)
	      {
			Q = min(0.9999999999, max(0.0000000001, Q));
		     return dH+n*(R*T*log((n-n*Q)*(1-Q)/((1+n*Q)*(n+Q))))/(n+1)+(2*Q-1)*W;
	      });
	
	double Q = 0.5*(a+b);

	double Xa1 = (1+n*Q)/(n+1);
	double Xs1 = (n-n*Q)/(n+1);
	double Xa2 = (1-Q)/(n+1);
	double Xs2 = (n+Q)/(n+1);

	double dS = -R*(Xa1*log(Xa1) + n*Xa2*log(Xa2) + Xs1*log(Xs1) + n*Xs2*log(Xs2));

	double Gdis = Fac*(dH + Q*(W - dH) - Q*Q*W - T*dS);
	
	Gfi += Gdis;
	
	return Gfi;
}

double Bragg_Williams_Solid::volume(Phase const &phase, double T, double P) const
{
  return phase.V;
}
