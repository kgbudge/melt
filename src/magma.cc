// magma.cc
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

/*#include <cmath>
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <vector>

#include "ds++/Assert.hh"
*/

#include "magma.hh"
#include "phase.hh"

//inline double square(double x){ return x*x; }

double Magma::Gf(Phase const &phase, double T, double P) const
{
	return phase.Hf0;
}

double Magma::volume(Phase const &phase, double T, double P) const
{
	return phase.V;
}
