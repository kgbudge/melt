/*
 * Melt_Model__mimimize_trial_set_.cc
 * Copyright (C) 2021 Kent G. Budge <kgb@kgbudge.com>
 * 
 * melt is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * melt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Melt_Model.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "algorithm.hh"

//-----------------------------------------------------------------------------//
double Melt_Model::minimize_trial_set_(unsigned const n,
                                       unsigned const cphase[P_END],
                                       double p[P_END])
{
	using namespace std;

	double g[n], h[n], xi[n];

	double xm[E_END];
	compute_current_melt_composition_(p, xm);

	unsigned const NPL = phase_.size();

	double fp = Gf(p);
	double pnorm = 0.0;
	for (unsigned i=0; i<n; ++i)
	{
		xi[i] = dGf(p, xm, fp, cphase[i]);
		// Prune to enforce constraints
		if (xi[i]>0.0)
		{
			if (p[cphase[i]]<=1.0e-9*cnorm_)
			{
				// Prune melt of solid phase not present
				xi[i] = 0.0; 
				p[cphase[i]] = 0.0;
			}
			else
			{
				cout << "Melting " << phase_[cphase[i]].name << endl;
			}
		}
		else if (xi[i]<0.0)
		{
			// Check for crystallization of element fully extracted
			Phase const &ph = phase_[cphase[i]];
			unsigned const Z = ph.nz;
			for (unsigned j=0; j<Z; ++j)
			{
				unsigned z = ph.z[j];
				if (xm[z]<1.0e-9*cnorm_)
				{
					xi[i] = 0.0;
				}
			}
			if (xi[i] != 0.0)
			{
				cout << "Crystallizing " << ph.name << endl;
			}
		}
		pnorm += xi[i]*xi[i];
	}
	pnorm = sqrt(pnorm);

	for (int j=0; j<n; ++j)
	{
		g[j] = -xi[j];
		xi[j] = h[j] = g[j];
	}

	for (;;)
	{
		if (pnorm < 2.0e-7*cnorm_) // To make sure search direction is meaningful
			return fp;

		// Compute freeze composition on pruned search direction. Clip upper
		// limit of search to not melt more of any phase than is actually present.
		// Also compute current melt composition.
		double xfreeze[E_END];
		fill(xfreeze, xfreeze+E_END, 0.0);
		double x0 = 0.0;
		double x1 = numeric_limits<double>::max();
		for (unsigned i=0; i<n; ++i)
		{
			unsigned const j = cphase[i];
			if (xi[i]!=0.0)
			{
				Phase const &ph = phase_[j];
				unsigned Z = ph.nz;
				for (unsigned k=0; k<Z; ++k)
				{
					unsigned z = ph.z[k];
					xfreeze[z] += xi[i]*ph.n[k];
				}
			}
			if (xi[i]<0.0)
			{
				x1 = min(x1, -p[j]/xi[i]);
			}
		}
		// Clip upper limit of search to not extract more freeze composition
		// than the amount of melt can supply.
		for (unsigned e=0; e<E_END; ++e)
		{
			if (xfreeze[e]>0.0)
			{
				x1 = min(x1, xm[e]/xfreeze[e]);
			}
		}

		double x2 = minimize(x0, x1, [&](double const e)
		                     {return Gfmelt(p,
		                                    n,
		                                    cphase,
		                                    xi,
		                                    e);});

		for (unsigned i=0; i<n; ++i)
		{
			unsigned const j = cphase[i];
			p[j] += x2*xi[i];
			p[j] = max(0.0, p[j]);
		}

		fp = this->Gf(p);
	    compute_current_melt_composition_(p, xm);

		cout << endl;
		cout << "Line search completed. New Gf = " << fp << endl;
		cout << "Solid composition after line search:" << endl;
		for (unsigned i=0; i<NPL; ++i)
		{
			if (p[i]>1.0e-9*cnorm_)
			{
				cout << "  " << phase_[i].name << " = " << p[i] << endl;
			}
		}
		cout << "Melt composition after line search:" << endl;
		for (unsigned i=0; i<E_END; ++i)
		{
			if (xm[i]>0.0)
			{
				cout << "  " << element_name[i] << " = " << xm[i] << endl;
			}
		}
		cout << endl;

		if (fabs(x2)*pnorm<3.0e-7*cnorm_) 
			return fp;

		cout << "New pruned gradient calculation:" << endl;
		double dgg = 0., gg = 0.;
		for (unsigned i=0; i<n; ++i)
		{
			unsigned const j = cphase[i];
			xi[i] = dGf(p, xm, fp, j);

			// Prune to enforce constraints
			if (xi[i]>0.0)
			{
				if (p[j]<=1.0e-9*cnorm_)
				{
					// Prune melt of solid phase not present
					p[j] = xi[i] = 0.0; 
					cout << phase_[j].name << " constrained from melting" << endl;
				}
				else
				{
					cout << "Melting " << phase_[j].name << ", p = " << p[j] << endl;
				}
			}
			else if (xi[i]<0.0)
			{
				// Check for crystallization of element fully extracted
				Phase const &ph = phase_[j];
				unsigned const Z = ph.nz;
				for (unsigned k=0; k<Z; ++k)
				{
					unsigned z = ph.z[k];
					if (xm[z]<1.0e-9*cnorm_)
					{
						xi[i] = 0.0;
				     	cout << ph.name << " constrained from crystallizing" << endl;
					}
				}
				if (xi[i] != 0.0)
				{
					cout << "Crystallizing " << ph.name << ", p = " << p[j]  << endl;
				}
			}
			gg += g[i]*g[i];
			dgg += (xi[i]+g[i])*xi[i];
		}

		if (gg == 0.0)
			break;

		double  gam = dgg/gg;
		pnorm = 0.0;
		cout << endl;
		cout << "Conjugate gradient search direction:" << endl;
		for (unsigned i=0; i<n; i++)
		{
			g[i] = -xi[i];
			if (g[i] != 0.0)
			{
		   	   xi[i] = h[i] = g[i] + gam*h[i];
			}
			else
			{
		   	   xi[i] = h[i] = 0.0;
			}

			if (xi[i]!= 0.0)
			{
				cout << phase_[cphase[i]].name << " p = " << xi[i] << endl;
				pnorm += xi[i]*xi[i];
			}
		}
		pnorm = sqrt(pnorm);
	}
	return fp;
}
