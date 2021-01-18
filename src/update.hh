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

#ifndef update_HH
#define update_HH

#include <string>
#include <vector>

#include "element.hh"
#include "phase.hh"
#include "State.hh"

void IUGS_classify(State const &);

bool 
do_ladder_update(double const T, 
            	 double const P, 
                 std::vector<Phase> const &phase,
                 std::vector<double> Gf,
                 bool const oxygen_specified, 
                 bool const oxygen_FMQ,
                 double &pO2,
                 State &);

void update_state(double const T, 
          		  double const P, 
           		  std::vector<Phase> &phase,
           		  bool const oxygen_specified, 
                  bool const oxygen_FMQ,
                  double &pO2,
	              State &state);

#endif // update_HH
