// Gf.cc
//
// Copyright (C) 2016 - Kent G. Budge
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
#include <algorithm>
#include <assert.h>

#include "ds++/Assert.hh"

#include "phase.hh"

inline double square(double x){ return x*x; }

using namespace std;

void initialize_Gf(double T /* K */, double P /* kbar */, double Gf[P_END] /* kJ */)
{
	// Note that Gf is the apparent free energy of formation, which ignores entropy of the elements.

	double const dT = T-T0;
	double const dT2 = T*T-T0*T0;
	double const drT = 1.0/T-1.0/T0;
	double const drT2 = 1.0/T/T - 1.0/T0/T0;
	double const sT = sqrt(T);
	double const sT0 = sqrt(T0);
	double const dsT = sT-sT0;
	double const drsT = 1.0/sT - 1.0/sT0;
	double const dlogT = log(T/T0);
	for (unsigned i=0; i<P_END; ++i)
	{
		double const Hf0 = phase[i].Hf0;
		double const S0 = phase[i].S0*1e-3; // to bring J to kJ
		double const V0 = phase[i].V;
		double const A = phase[i].a;
		double const B = phase[i].b*1e-5;  // By convention, reported in units of 1e-5 kJ/K/K
		double const C = phase[i].c;
		double const D = phase[i].d;
		double const a0 = phase[i].a0*1e-5; // By convention, reported in units of 1e-5/K
		double const k0 = phase[i].k0;
		double const k0p = phase[i].k0p;
		double const k0pp = phase[i].k0pp;
		unsigned const N = phase[i].nz;
		Model const model = phase[i].model;

		if (model == SOLID || model == MELT)
		{
			// Temperature terms
			double const Gt = Hf0 - T*S0 + A*dT + 0.5*B*dT2 - C*drT + 2*D*dsT
				- T*(A*dlogT + B*dT - 0.5*C*drT2 - 2*D*drsT);

			// Pressure term
			double const b = k0p*(2+k0p)/(k0*(1+k0p));
			double const c = 1.0/(k0p*(2+k0p));
			double const a = sqrt((1+c)/c);

			if (model == SOLID)
			{
				unsigned n = 0;
				for (unsigned j=0; j<phase[i].nz; ++j)
				{
					n += phase[i].n[j];
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

					Gf[i] = Gfi;
				}
				else
				{
					Gf[i] = 1.0e5; // off range of validity of table
				}
			}
			else  // model == MELT
			{
				double const Vt = V0*(1 + a0*dT - 20*a0*dsT);
					
				double const Gfi =
					Gt + 
					  P*Vt*(1-a+a*(1 - pow(1+b* P, 1-c))/(b*(c-1)* P))
					-P0*Vt*(1-a+a*(1 - pow(1+b*P0, 1-c))/(b*(c-1)*P0));


				Gf[i] = Gfi; 
			}
		}
		else if (model == VAPOR)
		{
			Gf[i] = solve_for_gibbs(phase[i], T, P);
		}
		else
		{
			Insist(false, "bad case");
		}
	}
}

double solve_for_fugacity(Phase const &g, double const T, double const Gf /* apparent free energy*/)
{
	double const a = g.a;
	double const b = g.b * 1.0e-5;
	double const c = g.c;
	double const d = g.d;
	double const S0 = g.S0 * 1e-3;
	double const Hf0 = g.Hf0;
	
	double const Gt = Hf0 - T*S0 + a*(T-T0) + 0.5*b*(T*T-T0*T0) - c*(1/T-1/T0) + 2*d*(sqrt(T)-sqrt(T0))
		- a*log(T/T0) + b*(T-T0) - 0.5*c*(1/(T*T)-1/(T0*T0)) - 2*d*(1/sqrt(T)-1/sqrt(T0)); 

	double const P = P0*exp((Gf - Gt)/(R*T));
	return P;
}

double solve_for_gibbs(Phase const &g, double const T, double const p)
{
	double const a = g.a;
	double const b = g.b * 1.0e-5;
	double const c = g.c;
	double const d = g.d;
	double const S0 = g.S0 * 1e-3;
	double const Hf0 = g.Hf0;
	
	double const Gt = Hf0 - T*S0 + a*(T-T0) + 0.5*b*(T*T-T0*T0) - c*(1/T-1/T0) + 2*d*(sqrt(T)-sqrt(T0))
		- a*log(T/T0) + b*(T-T0) - 0.5*c*(1/(T*T)-1/(T0*T0)) - 2*d*(1/sqrt(T)-1/sqrt(T0)); 

	double const Gf = Gt + R*T*log(p/P0);
	return Gf;
}
