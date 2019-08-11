// bm_solid.hh
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

#ifndef bm_solid_hh
#define bm_solid_hh

#include "solid.hh"

class BM_Solid : public Solid
{
  public:
    struct BM_Solid_Phase
    {
      double c0, dum1, c2, c3, c1; // heat capacity coeffs
      double dvdt, dvdp, d2vdpdt, Kp; // volume coeffs
    };

    virtual double Gf(Phase const &phase, double T, double P) const;
    virtual double volume(Phase const &phase, double T, double P) const;

  private:
    static double p(BM_Solid_Phase const &sph, double const Vl);

    static BM_Solid const *p_phase;	
    static double p_T, p_V;
};

#endif // bm_Solid_hh
