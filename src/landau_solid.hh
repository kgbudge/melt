// landau_solid.hh
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

#ifndef landau_solid_hh
#define landau_solid_hh

#include "solid.hh"

class Landau_Solid : public Solid
{
  struct Solid_Phase
  {
    double a, b, dm0, c, d; // heat capacity coeffs
    double a0, k0, k0p, k0pp; // volume coeffs
    double Tc, Smax, Vmax; // Landau coeffs
  };

  public:
    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;
};

#endif // landau_solid_hh
