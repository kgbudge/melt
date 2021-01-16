/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cc
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

#include "kgb/tostring.h"

#include "composition.hh"
#include "gui.hh"

//-----------------------------------------------------------------------------//
void on_normalize()
{
	double total = SiO2 + TiO2 + Al2O3 + Fe2O3 + FeO + MnO + MgO + CaO + Na2O + K2O + P2O5 + S + Cr2O3 + ZrO2 + H2O + CO2 + Cl;
	if (total>0.0)
	{
		double rnorm = 100.0/total;
		entry_SiO2->set_text(tostring(SiO2*rnorm));
		entry_TiO2->set_text(tostring(TiO2*rnorm));
		entry_Al2O3->set_text(tostring(Al2O3*rnorm));
		entry_Fe2O3->set_text(tostring(Fe2O3*rnorm));
		entry_FeO->set_text(tostring(FeO*rnorm));
		entry_MnO->set_text(tostring(MnO*rnorm));
		entry_MgO->set_text(tostring(MgO*rnorm));
		entry_CaO->set_text(tostring(CaO*rnorm));
		entry_Na2O->set_text(tostring(Na2O*rnorm));
		entry_K2O->set_text(tostring(K2O*rnorm));
		entry_P2O5->set_text(tostring(P2O5*rnorm));
		entry_S->set_text(tostring(S*rnorm));
		entry_Cr2O3->set_text(tostring(Cr2O3*rnorm));
		entry_ZrO2->set_text(tostring(ZrO2*rnorm));
		entry_H2O->set_text(tostring(H2O*rnorm));
		entry_CO2->set_text(tostring(CO2*rnorm));
		entry_Cl->set_text(tostring(Cl*rnorm));
		update();
	}
}
