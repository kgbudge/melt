// phase.hh
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

#ifndef phase_hh
#define phase_hh

#include "phase_enum.hh"

class Model;

unsigned const MAX_Z = 9;

extern struct Phase
{
    unsigned index; // should match enumerator
	char const *name;
	unsigned nz; // elements in formula
	unsigned z[MAX_Z]; // elements of formula
	double n[MAX_Z];  // quantities of each element in formula. double because these can be fractional for mineraloids or solid solutions.

	double Hf0;  // standard Gibbs free energy of formation at STP in kJ
	double S0; // entropy at STP in J/K
	double V; // molar volume at STP in kJ/kbar = 12.342 cm^3

	Model const *model;

	double data[60];
}
const phase[P_END];

#endif // phase_hh
