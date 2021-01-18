/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * classify.hh
 * Copyright (C) 2021 Kent G. Budge <kgb@kgbudge.com>
 * 
 * melt is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * melt is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef melt_classify_hh
#define melt_classify_hh

#include "State.hh"

//-----------------------------------------------------------------------------//
void classify_components(State const &state, 
double &Q, double &an, double &anc, double &lc, double &ne, double &An, 
double &A, double &P, double &F, 
double &M, double &Ol, double &Opx, double &Cpx,
std::string &m1, std::string &m2, std::string &um);

#endif

