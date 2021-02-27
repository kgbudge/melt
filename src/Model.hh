// model.hh
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

#ifndef model_hh
#define model_hh

struct Phase;

class Model
{
  // Temperature in K, pressure in kbar, Gibbs free energy in kJ/mol

  // Note that Gf is the apparent free energy of formation, which ignores the
  // entropy of the elements in their standard state.

  public:
    virtual double Gf(Phase const &phase, double T, double P) const = 0;
    virtual double volume(Phase const &phase, double T, double P) const = 0;
};

extern Model const *const SOLID;
extern Model const *const LANDAU_SOLID;
extern Model const *const BRAGG_WILLIAMS_SOLID;
extern Model const *const MELT;
extern Model const *const BM_MELT;
extern Model const *const BM_SOLID;
extern Model const *const VAPOR;
extern Model const *const PFITZER_STERNER;

#endif // model_hh
