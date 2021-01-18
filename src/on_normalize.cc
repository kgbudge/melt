/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * on_normalize.cc
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
void on_normalize()
{
	using namespace std;

	string name;
	bool is_molar;

	double SiO2;
	double TiO2;
	double Al2O3;
	double Fe2O3;
	double FeO;
	double MnO;
	double MgO;
	double CaO;
	double Na2O;
	double K2O;
	double P2O5;
	double S;
	double Cr2O3;
	double ZrO2;
	double H2O;
	double CO2;
	double Cl;

	get_composition(name,
	                is_molar,
	                SiO2,
	                TiO2,
	                Al2O3,
	                Fe2O3,
	                FeO,
	                MnO,
	                MgO,
	                CaO,
	                Na2O,
	                K2O,
	                P2O5,
	                S,
	                Cr2O3,
	                ZrO2,
	                H2O,
	                CO2,
	                Cl);

	// Mormalize to 100

	double const total = SiO2 + TiO2 + Al2O3 + Fe2O3 + FeO + MnO + MgO + 
		CaO + Na2O + K2O + P2O5 + S + Cr2O3 + ZrO2 + H2O + CO2 + Cl;

	if (total>0.0)
	{
		double rnorm = 100.0/total;

		entry_SiO2->set_text(tostring(rnorm*SiO2));
		entry_TiO2->set_text(tostring(rnorm*TiO2));
		entry_Al2O3->set_text(tostring(rnorm*Al2O3));
		entry_Fe2O3->set_text(tostring(rnorm*Fe2O3));
		entry_FeO->set_text(tostring(rnorm*FeO));
		entry_MnO->set_text(tostring(rnorm*MnO));
		entry_MgO->set_text(tostring(rnorm*MgO));
		entry_CaO->set_text(tostring(rnorm*CaO));
		entry_Na2O->set_text(tostring(rnorm*Na2O));
		entry_K2O->set_text(tostring(rnorm*K2O));
		entry_P2O5->set_text(tostring(rnorm*P2O5));
		entry_S->set_text(tostring(rnorm*S));
		entry_Cr2O3->set_text(tostring(rnorm*Cr2O3));
		entry_ZrO2->set_text(tostring(rnorm*ZrO2));
		entry_H2O->set_text(tostring(rnorm*H2O));
		entry_CO2->set_text(tostring(rnorm*CO2));
		entry_Cl->set_text(tostring(rnorm*Cl));
	}
}
