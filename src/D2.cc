// D2.cc
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

#include "D2.hh"

#include <cmath>

void D2::swap(D2 &r)
{
	std::swap(f,r.f);
	df.swap(r.df);
	ddf.swap(r.ddf);
}

D2 operator*(D2 const &a, D2 const &b)
{
	unsigned const N = a.df.size();
	D2 y(N);
	y.f = a.f*b.f;
	for (unsigned i=0; i<N; ++i)
	{
		Check(i<y.df.size() && i<b.df.size() && i<a.df.size());
		y.df[i] = a.f*b.df[i] + a.df[i]*b.f;

		Check(i<y.ddf.size() && i<b.ddf.size() && i<a.ddf.size());
		for (unsigned j=0; j<N; ++j)
		{
		    Check(j<y.ddf[i].size() && j<b.ddf[i].size() && j<a.ddf[i].size());
			y.ddf[i][j] = a.df[j]*b.df[i] + a.f*b.ddf[i][j] + a.ddf[i][j]*b.f + a.df[i]*b.df[j];
		}
	}
	return y;
}

D2 operator*(D2 const &a, double const b)
{
	D2 y(a);
	y *= b;
	return y;
}

D2 operator*(D2 &&a, double const b)
{
	D2 y;
	y.swap(a);
	y *= b;
	return y;
}

D2 &D2::operator*=(D2 const &r)
{
	// Safest implementation if not most efficient
	D2 y = *this * r;
	*this = y;
	return *this;
}

D2 &D2::operator*=(double r)
{
	unsigned N = df.size();
	f *= r;
	for (unsigned i=0; i<N; ++i)
	{
		Check(i<df.size());
		df[i] *= r;

		Check(i<ddf.size());
		for (unsigned j=0; j<N; ++j)
		{
			Check(j<ddf.size());
			ddf[i][j] *= r;
		}
	}
	return *this;
}

D2 operator/(D2 const &a, D2 const &b)
{
	unsigned N = b.df.size();
	D2 y(N);
	y.f = a.f/b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = a.df[i]/b.f - a.f*b.df[i]/(b.f*b.f);
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = a.ddf[i][j]/b.f - a.df[i]*b.df[j]/(b.f*b.f) - a.df[j]*b.df[i]/(b.f*b.f)
				- a.f*b.ddf[i][j]/(b.f*b.f) + 2*a.f*b.df[i]*b.df[j]/(b.f*b.f*b.f);
		}
	}
	return y;
}

D2 &D2::operator/=(D2 const &b)
{
	unsigned N = b.df.size();
	for (unsigned i=0; i<N; ++i)
	{
		for (unsigned j=0; j<N; ++j)
		{
			ddf[i][j] = ddf[i][j]/b.f - df[i]*b.df[j]/(b.f*b.f) - df[j]*b.df[i]/(b.f*b.f)
				- f*b.ddf[i][j]/(b.f*b.f) + 2*f*b.df[i]*b.df[j]/(b.f*b.f*b.f);
		}
	}
	for (unsigned i=0; i<N; ++i)
	{
		df[i] = df[i]/b.f - f*b.df[i]/(b.f*b.f);
	}
	f /= b.f;
	
	return *this;
}

D2 operator/(D2 &&a, D2 const &b)
{
	D2 y;
	y.swap(a);
    y /= b;
	return y;
}

D2 operator/(double const a, D2 const &b)
{
	unsigned N = b.df.size();
	D2 y(N);
	y.f = a/b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -a*b.df[i]/(b.f*b.f);
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = 2*a*b.df[i]*b.df[j]/(b.f*b.f*b.f) - a*b.ddf[i][j]/(b.f*b.f);
		}
	}
	return y;
}

D2 operator/(double const a, D2 &&b)
{
	unsigned N = b.df.size();
	D2 y;
	y.swap(b);
	for (unsigned i=0; i<N; ++i)
	{
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = 2*a*y.df[i]*y.df[j]/(y.f*y.f*y.f) - a*y.ddf[i][j]/(y.f*y.f);
		}
	}
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -a*y.df[i]/(y.f*y.f);
	}
	y.f = a/y.f;
	return y;
}

D2 operator+(D2 const &a, D2 const &b)
{
	unsigned const N = a.df.size();
	D2 y(N);
	y.f = a.f + b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = a.df[i] + b.df[i];
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = a.ddf[i][j] + b.ddf[i][j];
		}
	}
	return y;
}

D2 &D2::operator+=(D2 const &r)
{
	unsigned const N = df.size();
	f += r.f;
	for (unsigned i=0; i<N; i++)
	{
		df[i] += r.df[i];
		for (unsigned j=0; j<N; ++j)
		{
			ddf[i][j] += r.ddf[i][j];
		}
	}
	return *this;
}

D2 operator+(D2 const &a, double const b)
{
	D2 Result = a;
	Result.f += b;
	return Result;
}

D2 operator-(D2 const &a, D2 const &b)
{
	unsigned const N = a.df.size();
	D2 y(N);
	y.f = a.f - b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = a.df[i] - b.df[i];
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = a.ddf[i][j] - b.ddf[i][j];
		}
	}
	return y;
}

D2 operator-(D2 &&a, D2 const &b)
{
	D2 y;
	y.swap(a);
	y -= b;
	return y;
}

D2 operator-(double const a, D2 const &b)
{
	unsigned const N = b.df.size();
    D2 y(N);
	y.f = a - b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -b.df[i];
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = - b.ddf[i][j];
		}
	} 
	return y;
}


D2 &D2::operator-=(D2 const &r)
{
	unsigned const N = df.size();
	f -= r.f;
	for (unsigned i=0; i<N; i++)
	{
		df[i] -= r.df[i];
		for (unsigned j=0; j<N; ++j)
		{
			ddf[i][j] -= r.ddf[i][j];
		}
	}
	return *this;
}

D2 operator-(D2 const &x)
{
	unsigned const N = x.df.size();
	D2 y(N);
	y.f = -x.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -x.df[i];
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = -x.ddf[i][j];
		}
	}
	return y;
}

D2 log(D2 const &x)
{
	unsigned const N = x.df.size();
	D2 y(N);
	y.f = ::log(x.f);
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = x.df[i]/x.f;
		for (unsigned j=0; j<N; ++j)
		{
			y.ddf[i][j] = -x.df[i]*x.df[j]/(x.f*x.f) + x.ddf[i][j]/x.f;
		}
	}
	return y;
}

