/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * on_activate_depth.cc
 * Copyright (C) 2015 Kent G. Budge <kgb@kgbudge.com>
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

#include "kgb/tostring.h"

#include "gui.hh"

//-----------------------------------------------------------------------------//
void on_activate_depth()
{
	std::string text = entry_depth->get_text();
	double depth = atof(text.c_str());
	text = entry_over_rho->get_text();
	double rho = atof(text.c_str());
	double P = 1e5*1e-9*depth*rho*980.665 + 1.0e-3;
	entry_P->set_text(tostring(P, 3));
	
	update();
}
