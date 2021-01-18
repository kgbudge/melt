/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * get_composition.cc
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


#include "kgb/tostring.h"

#include "gui.hh"

//-----------------------------------------------------------------------------//
double 
get_composition(std::string &name,
                bool &is_molar,
	double &SiO2,
	double &TiO2,
	double &Al2O3,
	double &Fe2O3,
	double &FeO,
	double &MnO,
	double &MgO,
	double &CaO,
	double &Na2O,
	double &K2O,
	double &P2O5,
	double &S,
	double &Cr2O3,
	double &ZrO2,
	double &H2O,
	double &CO2,
	double &Cl)
{
	using namespace std;
	
	name = entry_name->get_text();
	is_molar = !button_byweight->get_active();

	string text = entry_SiO2->get_text();
	SiO2 = atof(text.c_str());
	text = entry_TiO2->get_text();
	TiO2 = atof(text.c_str());
	text = entry_Al2O3->get_text();
	Al2O3 = atof(text.c_str());
	text = entry_Fe2O3->get_text();
	Fe2O3 = atof(text.c_str());
	text = entry_FeO->get_text();
	FeO = atof(text.c_str());
	text = entry_MnO->get_text();
	MnO = atof(text.c_str());
	text = entry_MgO->get_text();
	MgO = atof(text.c_str());
	text = entry_CaO->get_text();
	CaO = atof(text.c_str());
	text = entry_Na2O->get_text();
	Na2O = atof(text.c_str());
	text = entry_K2O->get_text();
	K2O = atof(text.c_str());
	text = entry_P2O5->get_text();
	P2O5 = atof(text.c_str());
	text = entry_S->get_text();
	S = atof(text.c_str());
	text = entry_Cr2O3->get_text();
	Cr2O3 = atof(text.c_str());
	text = entry_ZrO2->get_text();
	ZrO2 = atof(text.c_str());
	text = entry_H2O->get_text();
	H2O = atof(text.c_str());
	text = entry_CO2->get_text();
	CO2 = atof(text.c_str());
	text = entry_Cl->get_text();
	Cl = atof(text.c_str());
	text = entry_T->get_text();
	double T = atof(text.c_str()) + 273.15;
	text = entry_P->get_text();
	double P = atof(text.c_str());

	double total = SiO2 + TiO2 + Al2O3 + Fe2O3 + FeO + MnO + MgO + CaO + Na2O +
		K2O + P2O5 + S + Cr2O3 + ZrO2 + H2O + CO2 + Cl;
	
	text_total->set_text(tostring(total));

	return total;
}
