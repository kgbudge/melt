/*
 * algorithm.hh
 * Copyright (C) 2019 Kent G. Budge <kgb@kgbudge.com>
 * 
 * norm is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * norm is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef algorithm_hh
#define algorithm_hh

#include <algorithm>
/*#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"

#include "element.hh"
#include "model.hh"
#include "phase.hh"*/

//-----------------------------------------------------------------------------//
using namespace std;

template<class F>
double minimize(double a, double b, F const &f)
{
	double fa = f(a), fb = f(b);
	double const phi = 0.5*(1+sqrt(5.));

	double c = b - (b-a)/phi;
	double fc =  f(c);
	double d = a + (b-a)/phi;
	double fd = f(d);

	while (fabs(b-a)>2.0e-12*(fabs(a)+fabs(b)))
	{
		if (fc<fd)
		{
			b = d;
			fb = fd;
			d = c;
			fd = fc;
			c = b - (b-a)/phi;
			fc =  f(c);
		}
		else
		{
			a = c;
			fa = fc;
			c = d;
			fc = fd;
			d = a + (b-a)/phi;
			fd = f(d);
		}
	};
	return (fb<fa? b : a);
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
	while (fb != 0 && fs != 0 && fabs(b-a)>1.0e-14*max(fabs(a),fabs(b)));
	x = fabs(fb)<=fabs(fs)? b : s;	
}

template<class Function>
void solvebr(double &a, double &b, Function f)
{
	// Root already bracketed.
	double fa = f(a), fb = f(b);
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

	    fs = f(s);
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
}

#endif // algorithm_hh
