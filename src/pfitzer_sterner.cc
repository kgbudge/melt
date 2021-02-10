// pfitzer_sterner.cc
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
#include <algorithm>
//#include <assert.h>
//#include <fstream>
//#include <vector>

#include "ds++/Assert.hh"

#include "algorithm.hh"
#include "constants.hh"
#include "pfitzer_sterner.hh"
#include "Phase.hh"

using namespace std;

inline double square(double x){ return x*x; }

static Pfitzer_Sterner const myPfitzer_Sterner;
Model const *const PFITZER_STERNER = &myPfitzer_Sterner;

double Pfitzer_Sterner::C[10];
double Pfitzer_Sterner::T;

double Pfitzer_Sterner::Gf(Phase const &phase, double const T, double const P) const
{
	// Calculate temperature-dependent coefficients in the Pfitzer-Stern
	// nonideal pressure formula
	
	for (unsigned i=0; i<10; ++i)
	{
	  C[i] = phase.data[i*6+3] + T*(phase.data[i*6+4] + T*phase.data[i*6+5]) 
			+ (phase.data[i*6+2] + (phase.data[i*6+1] + phase.data[i*6+0]/(T*T))/T)/T;
	}
    Pfitzer_Sterner::T = T;
	
	// Calculate nonideal molar density in mol/cm^3
    double rho = Pfitzer_Sterner::rho(P);

	// Calculate the ideal Gibbs free energy for this temperature and density.

	double G_ideal = Vapor::Gf(phase, T, P);
		
	double const Art_resid = C[0]*rho + (1/(C[1]+rho*(C[2]+rho*(C[3]+rho*(C[4]+rho*C[5])))) - 1/C[1])
	                   -(C[6]/C[7])*expm1(-C[7]*rho)-(C[8]/C[9])*expm1(-C[9]*rho);
		
	double const A_resid = R*T*Art_resid;
	
    return G_ideal + A_resid;
}

double Pfitzer_Sterner::volume(Phase const &phase, double const T, double const P) const
{
	return  0.1/Pfitzer_Sterner::rho(phase, T, P);   // rho M/ml -> M/dl
}

double Pfitzer_Sterner::p(double const rho)
{
	double const Prt = rho*(1 + rho*(C[0] - ((C[2]+rho*(2*C[3] + rho*(3*C[4] + rho*4*C[5])))/square(C[1]+rho*(C[2]+rho*(C[3]+rho*(C[4]+rho*C[5])))))
	                               + C[6]*exp(-C[7]*rho) + C[8]*exp(-C[9]*rho)));
	
	return Prt*R*T*10;  // rho -> mol/cc to mol/dl
}

double Pfitzer_Sterner::rho(Phase const &phase, double T, double P) const
{

	for (unsigned i=0; i<10; ++i)
	{
	  C[i] = phase.data[i*6+3] + T*(phase.data[i*6+4] + T*phase.data[i*6+5]) 
			+ (phase.data[i*6+2] + (phase.data[i*6+1] + phase.data[i*6+0]/(T*T))/T)/T;
	}
	Pfitzer_Sterner::T = T;

	return rho(P);
}

double Pfitzer_Sterner::rho(double P)
{
	double rho = min(1.0/18.0, P/(R*T*10));  // from mol/dl to mol/cc
    solve(P, rho, Pfitzer_Sterner::p);
	return rho;
}
