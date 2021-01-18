/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * do_save.cc
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

#include <fstream>

#include "gui.hh"

//-----------------------------------------------------------------------------//
void do_save()
{
	using namespace std;

	ofstream out(file.c_str());
	if (!out)
	{
		Gtk::MessageDialog dialog(*main_win, "Save Fail");
		dialog.set_secondary_text("Could not write save file " + file);

		dialog.run();
	}
	else
	{
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

		out << "name " << name << endl;
		if (is_molar)
		{
			out << "molar" << endl;
		}
		out << "SiO2 " << SiO2 << endl;
		out << "TiO2 " << TiO2 << endl;
		out << "Al2O3 " << Al2O3 << endl;
		out << "Fe2O3 " << Fe2O3 << endl;
		out << "FeO " << FeO << endl;
		out << "MnO " << MnO << endl;
		out << "MgO " << MgO << endl;
		out << "CaO " << CaO << endl;
		out << "Na2O " << Na2O << endl;
		out << "K2O " << K2O << endl;
		out << "P2O5 " << P2O5 << endl;
		out << "S " << S << endl;
		out << "Cr2O3 " << Cr2O3 << endl;
		out << "ZrO2 " << ZrO2 << endl;
		out << "H2O " << H2O << endl;
		out << "CO2 " << CO2 << endl;
		out << "Cl " << Cl << endl;
	}
}
