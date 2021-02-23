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
                                       Reaction const cphase[P_END],
                                       double p[P_END])
{
	using namespace std;

	double g[n], h[n], xi[n];
	double xm[E_END];

	compute_current_melt_composition_(p, xm);

	double fp = Gf(p);
	double fm = Gfm(xm);
	double pnorm = 0.0;
	cout << "Active reactions:" << endl;
	for (unsigned i=0; i<n; ++i)
	{
		xi[i] = dGfm(p, xm, fm, cphase[i].i) + cphase[i].dGfs;
		// Prune to enforce constraints
		double r = calculate_extent_(cphase[i], p, xm);
		if (xi[i]>0.0)
		{
			if (r<=1.0e-9*cnorm_)
			{
				// Prune melt of solid phase not present
				xi[i] = 0.0; 
			}
			else
			{
				print_reverse(cphase[i]);
			}
		}
		else if (xi[i]<0.0)
		{
			if (r<=1.0e-9*cnorm_)
			{
				// Prune crystallization of solid phase not present
				xi[i] = 0.0; 
			}
			else
			{
				print(cphase[i]);
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

		// Determine change of solid and melt compositions per unit of search.
		double xmelt[E_END];
		fill(xmelt, xmelt+E_END, 0.0);
		double xsolid[P_END];
		fill(xsolid, xsolid+NP_, 0.0);
		for (unsigned i=0; i<n; ++i)
		{
			if (xi[i]!=0.0)
			{
				Reaction const &react = cphase[i];
				Phase const &ph = phase_[react.i];
				unsigned Z = ph.nz;
				for (unsigned k=0; k<Z; ++k)
				{
					xmelt[ph.z[k]] -= xi[i]*ph.n[k];
				}
				unsigned const N = react.nz;
				for (unsigned j=0; j<N; ++j)
				{
					xsolid[react.p[j]] += xi[i]*react.n[j];
				}
			}
		}
		// Clip upper limit of search to not extract more melt
		// or any solid phase than the current composition permits.
		double x0 = 0.0;
		double x1 = numeric_limits<double>::max();
		for (unsigned e=0; e<E_END; ++e)
		{
			if (xmelt[e]<0.0)
			{
				x1 = min(x1, -xm[e]/xmelt[e]);
			}
		}
		for (unsigned e=0; e<NP_; ++e)
		{
			if (xsolid[e]<0.0)
			{
				x1 = min(x1, -p[e]/xsolid[e]);
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
			auto const &r = cphase[i];
			unsigned const N = r.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned const z = r.p[j];
		  	    p[z] += x2*xi[i]*r.n[j];
			    p[z] = max(0.0, p[z]);
			}
		}

		fp = this->Gf(p);
	    compute_current_melt_composition_(p, xm);
		fm = Gfm(xm);

		cout << endl;
		cout << "Line search completed. New Gf = " << fp << endl;
		cout << "Solid composition after line search:" << endl;
		for (unsigned i=0; i<NP_; ++i)
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
			auto const r = cphase[i];
			xi[i] = dGfm(p, xm, fm, r.i) + r.dGfs;

			// Prune to enforce constraints
			if (xi[i]>0.0)
			{
				/*
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
			*/
			}
			else if (xi[i]<0.0)
			{
				// Check for crystallization of element fully extracted
				double ext = calculate_extent_(r, p, xm);
				if (ext<1.0e-9*cnorm_)
				{
					cout << phase_[r.i].name << " constrained from crystallizing" << endl;
					xi[i] = 0.0;
				}
				else
				{
					cout << "Crystallizing " << phase_[r.i].name << endl;
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
//				cout << phase_[cphase[i]].name << " p = " << xi[i] << endl;
				pnorm += xi[i]*xi[i];
			}
		}
		pnorm = sqrt(pnorm);
	}
	return fp;
}
