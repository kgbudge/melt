/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * on_activate_P.cc
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
void on_activate_P()
{
	std::string text = entry_P->get_text();
	double P = atof(text.c_str());
	text = entry_over_rho->get_text();
	double rho = atof(text.c_str());
	double depth = 1e-5*1e9*(P-1e-3)/(rho*980.665);
	entry_depth->set_text(tostring(depth, 3));
	
	update();
}
