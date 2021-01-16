/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * on_new.cc
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

#include "gui.hh"

//-----------------------------------------------------------------------------//
void on_new()
{
	file = "";
	name = "";
	
	entry_name->set_text("");
	entry_SiO2->set_text("0.00");
	entry_TiO2->set_text("0.00");
	entry_Al2O3->set_text("0.00");
	entry_Fe2O3->set_text("0.00");
	entry_FeO->set_text("0.00");
	entry_MnO->set_text("0.00");
	entry_MgO->set_text("0.00");
	entry_CaO->set_text("0.00");
	entry_Na2O->set_text("0.00");
	entry_K2O->set_text("0.00");
	entry_P2O5->set_text("0.00");
	entry_S->set_text("0.00");
	entry_Cr2O3->set_text("0.00");
	entry_ZrO2->set_text("0.00");
	entry_H2O->set_text("0.00");
	entry_CO2->set_text("0.00");
	entry_Cl->set_text("0.00");
	
	update();
}
