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
	for (unsigned i=0; i<P_END; ++i)
	{
      Gf[i] = phase[i].model->Gf(phase[i], T, P);
	}

	// Water critical cutoff
	if (T>674.096)
	{
		  Gf[P_H2O_LIQUID] = 1e5;
	}
}

//------------------------- Model ----------------------------------------------

double Model::P(Phase const &g, double const T, double const Gf) const
{
	  Insist(false, "not implemented");
}

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

	double const Hf0 = phase.Hf0;
	double const S0 = phase.S0*1e-3; // to bring J to kJ
	double const V0 = phase.V;
	double const A = phase.a;
	double const B = phase.b*1e-5;  // By convention, reported in units of 1e-5 kJ/K/K
	double const C = phase.c;
	double const D = phase.d;
	double const a0 = phase.a0*1e-5; // By convention, reported in units of 1e-5/K
	double const k0 = phase.k0;
	double const k0p = phase.k0p;
	double const k0pp = phase.k0pp;
	unsigned const N = phase.nz;

	// Temperature terms
	double const Gt = Hf0 - T*S0 + A*dT + 0.5*B*dT2 - C*drT + 2*D*dsT
		- T*(A*dlogT + B*dT - 0.5*C*drT2 - 2*D*drsT);

	// Pressure term
	double const b = k0p*(2+k0p)/(k0*(1+k0p));
	double const c = 1.0/(k0p*(2+k0p));
	double const a = sqrt((1+c)/c);

	unsigned n = 0;
	for (unsigned j=0; j<phase.nz; ++j)
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

//------------------------- Vapor ----------------------------------------------

double Vapor::Gf(Phase const &phase, double const T, double const P) const
{
	// Note that Gf is the apparent free energy of formation (in kJ/mol),
	// which ignores entropy of the elements.

	// Temperature is in K and P in kbar.

	double const a = phase.a;
	double const b = phase.b * 1.0e-5;
	double const c = phase.c;
	double const d = phase.d;
	double const S0 = phase.S0 * 1e-3;
	double const Hf0 = phase.Hf0;

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
	double const a = phase.a;
	double const b = phase.b * 1.0e-5;
	double const c = phase.c;
	double const d = phase.d;
	double const S0 = phase.S0 * 1e-3;
	double const Hf0 = phase.Hf0;
	
	double const Gt = Hf0 - T*S0 + a*(T-T0) + 0.5*b*(T*T-T0*T0) - c*(1/T-1/T0) + 2*d*(sqrt(T)-sqrt(T0))
		- a*log(T/T0) + b*(T-T0) - 0.5*c*(1/(T*T)-1/(T0*T0)) - 2*d*(1/sqrt(T)-1/sqrt(T0)); 

	double const P = P0*exp((Gf - Gt)/(R*T));
	return P;
}

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

	double const Hf0 = phase.Hf0;
	double const S0 = phase.S0*1e-3; // to bring J to kJ
	double const V0 = phase.V;
	double const A = phase.a;
	double const B = phase.b*1e-5;  // By convention, reported in units of 1e-5 kJ/K/K
	double const C = phase.c;
	double const D = phase.d;
	double const a0 = phase.a0*1e-5; // By convention, reported in units of 1e-5/K
	double const k0 = phase.k0;
	double const k0p = phase.k0p;
	double const k0pp = phase.k0pp;
	unsigned const N = phase.nz;

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

//------------------------- Aqueous ---------------------------------------------

double Aqueous::Gf(Phase const &phase, double const T, double const P) const
{
	Insist(false, "not yet implemented");
}

double Aqueous::volume(Phase const &phase, double T, double P) const
{
	Insist(false, "not yet implemented");
}


//------------------------- Water ----------------------------------------------
double Water::c1, Water::c2, Water::c3, Water::c4, Water::c5, Water::c6;
double Water::c7, Water::c8, Water::c9, Water::c0, Water::T;

double Water::Gf(Phase const &phase, double const T, double const P) const
{
	// Calculate temperature-dependent coefficients in the Pfitzer-Stern
	// nonideal pressure formula.
	c1 = c13/T + c14;   
	c2 = c23/T + c24 + c25*T;
	c3 = c33/T + c34 + T*(c35 + T*c36);
	c4 = c44 + T*c45;
	c5 = c53/T + c54 + T*c55;
	c6 = c64;
	c7 = (c73 + (c72 + c71/(T*T))/T)/T + c74;
	c8 = c83/T + c84;
	c9 = (c93 + (c92 + c91/(T*T))/T)/T + c94;
	c0 = c03/T + c04;
	Water::T = T;
	
	// Calculate nonideal molar density in mol/cm^3
    double rho = Water::rho(P);

	// Calculate the ideal Gibbs free energy for this temperature and density.

	double G_ideal = Vapor::Gf(phase, T, P);
		
	double const Art_resid = c1*rho + (1/(c2+rho*(c3+rho*(c4+rho*(c5+rho*c6)))) - 1/c2)
	                   -(c7/c8)*expm1(-c8*rho)-(c9/c0)*expm1(-c0*rho);
		
	double const A_resid = R*T*Art_resid;
	
    return G_ideal + A_resid;
}

double Water::volume(Phase const &phase, double const T, double const P) const
{
	return  0.1/Water::rho(phase, T, P);   // rho M/ml -> M/dl
}

double Water::p(double const rho)
{
	double const Prt = rho*(1 + rho*(c1 - ((c3+rho*(2*c4 + rho*(3*c5 + rho*4*c6)))/square(c2+rho*(c3+rho*(c4+rho*(c5+rho*c6)))))
	                               + c7*exp(-c8*rho) + c9*exp(-c0*rho)));
	
	return Prt*R*T*10;  // rho -> mol/cc to mol/dl
}

double Water::rho(Phase const &phase, double T, double P)
{

	c1 = c13/T + c14;   
	c2 = c23/T + c24 + c25*T;
	c3 = c33/T + c34 + T*(c35 + T*c36);
	c4 = c44 + T*c45;
	c5 = c53/T + c54 + T*c55;
	c6 = c64;
	c7 = (c73 + (c72 + c71/(T*T))/T)/T + c74;
	c8 = c83/T + c84;
	c9 = (c93 + (c92 + c91/(T*T))/T)/T + c94;
	c0 = c03/T + c04;
	Water::T = T;

	return rho(P);
}

double Water::rho(double P)
{
	double rho = min(1.0/18.0, P/(R*T*10));  // from mol/dl to mol/cc
    solve(P, rho, Water::p);
	return rho;
}

template<class Function>
void solve(double const y, double &x, Function f)
{
	// Bracket root. We know it's a monotonic increasing function.
	double a = 0.98*x, b = 1.02*x;
	double fa = f(a)-y, fb = f(b)-y;
	while (fa*fb>=0.0)
	{
		if (fabs(fa)<fabs(fb))
		{
			a *= 0.5;
			fa = f(a)-y;
			b *= 1.05;
			fb = f(b)-y;
		}
		else
		{
			b *= 1.4;
			fb = f(b)-y;
			a *= 0.95;
			fa = f(a)-y;
		}
	}

	// Now that root is bracketed, close in on it
	    
    if (fabs(fa) < fabs(fb)) 
	{
      swap (a,b);
  	  swap(fa,fb);
	}
	double c = a;
	double fc = fa;
	bool mflag = true;
	double s, fs, d;
	do {
		if (fa != fc && fb != fc)
		{
			s = (a*fb*fc)/((fa-fb)*(fa-fc))+(b*fa*fc)/((fb-fa)*(fb-fc))+(c*fa*fb)/((fc-fa)*(fc-fb));
		}
		else
		{
			s = b-fb*(b-a)/(fb-fa);
		}
		if ((s-(3*a+b)/4)*(s-b)>=0 ||
		    mflag && fabs(s-b) >= fabs(b-c)/2 ||
		    !mflag && fabs(s-b) >= fabs(c-d)/2 ||
		    mflag && fabs(b-c) < 1.0e-7*fabs(a-b) ||
		    !mflag && fabs(c-d) < 1.0e-7*fabs(a-b))
		{
			s = 0.5*(a+b);
			mflag = true;
		}
		else
		{
			mflag = false;
		}

	    fs = f(s) - y;
		d = c;
		c = b;
		fc = fb;
		if (fa*fs<0)
		{
			b = s;
			fb = fs;
		}
		else
		{
			a = s;
			fa = fs;
		}
		if (fabs(fa) < fabs(fb))
		{
			swap(a,b);
			swap(fa, fb);
		}
	}
	while (fb != 0 && fs != 0 && fabs(b-a)>1.0e-14*max(b,a));
	x = fabs(fb)<=fabs(fs)? b : s;	

#if 0
	double x2 = (y - y0)*(x1-x0)/(y1-y0) + x0;	
	while (x1>x0 && x2 != x1 && x2 != x0) 
	{
		double y2  = f(x2);
		if ((y2>y)==(y0>y))
		{
			y0 = y2;
			x0 = x2;
		}
		else if ((y2>y)==(y1>y))
		{
			y1 = y2;
			x1 = x2;
		}
		else
		{
			Insist(false, "construction");
		}
	    x2 = (y - y0)*(x1-x0)/(y1-y0) + x0;	
	}
	x = x2;
#endif
}

