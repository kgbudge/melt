// D1.cc
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

#include "D1.hh"

#include <cmath>

void D1::swap(D1 &r)
{
	std::swap(f,r.f);
	df.swap(r.df);
}

D1 operator*(D1 const &a, D1 const &b)
{
	unsigned const N = a.df.size();
	D1 y(N);
	y.f = a.f*b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = a.f*b.df[i] + a.df[i]*b.f;
	}
	return y;
}

D1 operator*(D1 const &a, double const b)
{
	D1 y(a);
	y *= b;
	return y;
}

D1 operator*(D1 &&a, double const b)
{
	D1 y;
	y.swap(a);
	y *= b;
	return y;
}

D1 &D1::operator*=(D1 const &r)
{
	// Safest implementation if not most efficient
	D1 y = *this * r;
	*this = y;
	return *this;
}

D1 &D1::operator*=(double r)
{
	unsigned N = df.size();
	f *= r;
	for (unsigned i=0; i<N; ++i)
	{
		Check(i<df.size());
		df[i] *= r;
	}
	return *this;
}

D1 operator/(D1 const &a, D1 const &b)
{
	unsigned N = b.df.size();
	D1 y(N);
	y.f = a.f/b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = a.df[i]/b.f - a.f*b.df[i]/(b.f*b.f);
	}
	return y;
}

D1 &D1::operator/=(D1 const &b)
{
	unsigned N = b.df.size();
	for (unsigned i=0; i<N; ++i)
	{
		df[i] = df[i]/b.f - f*b.df[i]/(b.f*b.f);
	}
	f /= b.f;
	
	return *this;
}

D1 operator/(D1 &&a, D1 const &b)
{
	D1 y;
	y.swap(a);
    y /= b;
	return y;
}

D1 operator/(double const a, D1 const &b)
{
	unsigned N = b.df.size();
	D1 y(N);
	y.f = a/b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -a*b.df[i]/(b.f*b.f);
	}
	return y;
}

D1 operator/(double const a, D1 &&b)
{
	unsigned N = b.df.size();
	D1 y;
	y.swap(b);
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -a*y.df[i]/(y.f*y.f);
	}
	y.f = a/y.f;
	return y;
}

D1 operator+(D1 const &a, D1 const &b)
{
	unsigned const N = a.df.size();
	D1 y(N);
	y.f = a.f + b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = a.df[i] + b.df[i];
	}
	return y;
}

D1 &D1::operator+=(double r)
{
	f += r;
	return *this;
}

D1 &D1::operator+=(D1 const &r)
{
	unsigned const N = df.size();
	f += r.f;
	for (unsigned i=0; i<N; i++)
	{
		df[i] += r.df[i];
	}
	return *this;
}

D1 operator+(D1 const &a, double const b)
{
	D1 Result = a;
	Result.f += b;
	return Result;
}

D1 operator-(D1 const &a, D1 const &b)
{
	unsigned const N = a.df.size();
	D1 y(N);
	y.f = a.f - b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = a.df[i] - b.df[i];
	}
	return y;
}

D1 operator-(D1 &&a, D1 const &b)
{
	D1 y;
	y.swap(a);
	y -= b;
	return y;
}

D1 operator-(double const a, D1 const &b)
{
	unsigned const N = b.df.size();
    D1 y(N);
	y.f = a - b.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -b.df[i];
	} 
	return y;
}


D1 &D1::operator-=(D1 const &r)
{
	unsigned const N = df.size();
	f -= r.f;
	for (unsigned i=0; i<N; i++)
	{
		df[i] -= r.df[i];
	}
	return *this;
}

D1 operator-(D1 const &x)
{
	unsigned const N = x.df.size();
	D1 y(N);
	y.f = -x.f;
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = -x.df[i];
	}
	return y;
}

D1 log(D1 const &x)
{
	unsigned const N = x.df.size();
	D1 y(N);
	y.f = ::log(x.f);
	for (unsigned i=0; i<N; ++i)
	{
		y.df[i] = x.df[i]/x.f;
	}
	return y;
}

