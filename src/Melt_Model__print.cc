/*
 * Melt_Model__print.cc
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

#include <iostream>

//-----------------------------------------------------------------------------//
void Melt_Model::print(Reaction const &r) const
{
	using namespace std;

	cout << phase_[r.i].name << ": melt ";
	cout << defaultfloat;
	unsigned N = r.nz;
	for (unsigned i=0; i<N; ++i)
	{
		if (r.n[i]<0.0)
		{
			cout << " + " << -r.n[i] << phase_[r.p[i]].name;
		}
	}
	cout << " -> ";
	bool first = true;
	for (unsigned i=0; i<N; ++i)
	{
		if (r.n[i]>0.0)
		{
			if (!first)
			{
				cout << " + ";
			}
			cout << r.n[i] << phase_[r.p[i]].name;
			first = false;
		}
	}
	cout << endl;
}

//-----------------------------------------------------------------------------//
void Melt_Model::print_reverse(Reaction const &r) const
{
	using namespace std;

	cout << phase_[r.i].name << ": ";
	cout << defaultfloat;
	unsigned N = r.nz;
	bool first = true;
	for (unsigned i=0; i<N; ++i)
	{
		if (r.n[i]>0.0)
		{
			if (!first)
			{
				cout << " + ";
			}
			cout << " + " << r.n[i] << phase_[r.p[i]].name;
			first = false;
		}
	}
	cout << " -> melt + ";
	for (unsigned i=0; i<N; ++i)
	{
		if (r.n[i]<0.0)
		{
			cout << " + ";
			cout << -r.n[i] << phase_[r.p[i]].name;
		}
	}
	cout << endl;
}
