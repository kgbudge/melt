/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * update.hh
 * Copyright (C) 2015 Kent G. Budge <kgb@kgbudge.com>
 * 
 * norm is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * norm is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef State_HH
#define State_HH

#include <string>

#include "element.hh"

// Deliberately keep this simple, so copy is relatively inexpensive in ladder
// The full description of the calculation requires an additional phase table
// and the T and P for which this calculation is done.
struct State
{
    std::string name;
	unsigned p[E_END];
	double x[E_END];
    bool is_element_active[E_END];
	double element_activity[E_END];
	double V[E_END];

	State(std::string const &name,        
          double nH2O, double nCO2, double nNa2O, double nMgO, double nAl2O3,
          double nSiO2, double nP2O5, double nS, double nCl, double nK2O, 
          double nCaO, double nTiO2, double nCr2O3, double nMnO, double nFeO, 
          double nFe2O3,  double nZrO2);

};

#endif // State_HH
