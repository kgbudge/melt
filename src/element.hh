// element.hh
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

#ifndef element_hh
#define element_hh

enum ELEMENT
{
  E_H,
  E_C,
  E_O,
  E_NA,
  E_MG,
  E_AL,
  E_SI,
  E_S,
  E_CL,
  E_K,
  E_CA,
  E_TI,
  E_CR,
  E_MN,
  E_FE,
  E_ZR,

  E_END
};

extern double const atomic_weight[E_END];

#endif // element_hh
