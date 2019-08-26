// element.cc
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

#include "element.hh"

extern double const atomic_weight[E_END] =
{
	1.008,
	12.011,
	15.999,
	22.990,
	24.305,
	26.982,
	28.085,
	30.974,
	32.06,
	35.45,
	39.098,
	40.078,
	47.867,
	51.996,
	54.938,
	55.845,
	91.224 
};

char const *const element_name[E_END] =
{
	"H",
	"C",
	"O",
	"Na",
	"Mg",
	"Al",
	"Si",
	"P",
	"S",
	"Cl",
	"K",
	"Ca",
	"Ti",
	"Cr",
	"Mn",
	"Fe",
	"Zr"
};

