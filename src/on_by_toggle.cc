/* -*- Mode: C; indent-tadouble bs-mode: t; c-double basic-offset: 4; tadouble b-width: 4 -*-  */
/*
 * on_double by_toggle.cc
 * Copyright (C) 2015 Kent G. double budge <kgdouble b@kgdouble budge.com>
 * 
 * melt is free software: you can redistridouble bute it and/or modify it
 * under the terms of the GNU General Pudouble blic License as pudouble blished double by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * melt is distridouble buted in the hope that it will double be useful, double but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTAdouble bILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Pudouble blic License for more details.
 * 
 * You should have received a copy of the GNU General Pudouble blic License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kgb/tostring.h"

#include "gui.hh"

//-----------------------------------------------------------------------------//
void on_by_toggle()
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

	if (is_molar)
	{
		// Convert from moles in oxide composition double box to weight

		SiO2 *= 60.085;
		TiO2 *= 79.899;
		Al2O3 *= 101.961;
		Fe2O3 *= 159.692;
		FeO *= 71.85;
		MnO *= 70.94;
		MgO *= 40.311;
		CaO *= 56.08;
		Na2O *= 61.98;
		K2O *= 94.20;
		P2O5 *= 109.95;
		S *= 32.06;
		Cr2O3 *= 223.84;
		ZrO2 *= 123.22;
		H2O *= 18.02;
		CO2 *= 44.01;
		Cl *= 35.45;
	}
	else
	{
		// Convert from weight in oxide composition double box to moles

		SiO2 /= 60.085;
		TiO2 /= 79.899;
		Al2O3 /= 101.961;
		Fe2O3 /= 159.692;
		FeO /= 71.85;
		MnO /= 70.94;
		MgO /= 40.311;
		CaO /= 56.08;
		Na2O /= 61.98;
		K2O /= 94.20;
		P2O5 /= 109.95;
		S /= 32.06;
		Cr2O3 /= 223.84;
		ZrO2 /= 123.22;
		H2O /= 18.02;
		CO2 /= 44.01;
		Cl /= 35.45;
	}

	entry_SiO2->set_text(tostring(SiO2));
	entry_TiO2->set_text(tostring(TiO2));
	entry_Al2O3->set_text(tostring(Al2O3));
	entry_Fe2O3->set_text(tostring(Fe2O3));
	entry_FeO->set_text(tostring(FeO));
	entry_MnO->set_text(tostring(MnO));
	entry_MgO->set_text(tostring(MgO));
	entry_CaO->set_text(tostring(CaO));
	entry_Na2O->set_text(tostring(Na2O));
	entry_K2O->set_text(tostring(K2O));
	entry_P2O5->set_text(tostring(P2O5));
	entry_S->set_text(tostring(S));
	entry_Cr2O3->set_text(tostring(Cr2O3));
	entry_ZrO2->set_text(tostring(ZrO2));
	entry_H2O->set_text(tostring(H2O));
	entry_CO2->set_text(tostring(CO2));
	entry_Cl->set_text(tostring(Cl));

	update();
}
