/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * classify_components.cc
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

#include <iostream>

#include "classify.hh"
#include "Phase.hh"

void sort_mafic(double x, std::string const &name,
                double &x1, std::string &m1,
                double &x2, std::string &m2)
{
	if (x>0.0)
	{
		if (x>x1)
		{
			m2.swap(m1);
			x2 = x1;
			x1 = x;
			m1 = name;
		}
		else if (x>x2)
		{
			x2 = x;
			m2 = name;
		}
	}
}

void sort_ultramafic(double x, std::string const &name,
                double &x1, std::string &m1)
{
	if (x>0.0)
	{
		if (x>x1)
		{
			x1 = x;
			m1 = name;
		}
	}
}

//-----------------------------------------------------------------------------//
void classify_components(State const &state, 
                         double &Q, double &an, double &anc, double &ne,
                         double &lc,
                         double &An, double &A, double &P, double &F, double &M,
                         double &ol, double &Opx, double &Cpx,
                         std::string &m1, std::string &m2, std::string &um)
{
	using namespace std;

	Q = an = anc = ne = lc = An = A = P = F = M = 0.0, ol = 0.0;
	Opx = Cpx = 0.0;
	double ab = 0.0, anda = 0.0, bio = 0, co = 0.0;
	double cor = 0, di = 0.0, ho = 0.0, hy = 0.0, mc = 0.0;
	double garnet = 0, spinel = 0;

	double const *const state_V = state.V();
	int const *const state_ph = state.ph();
	auto const &state_phase = state.phase();
	for (unsigned i=0; i<E_END; ++i)
	{
		double const x = state_V[i];
		int idx = state_phase[state_ph[i]].index;
		if (x>0.0 && idx<P_END)
		{
			switch (idx)
			{				
				case P_ACMITE:
					M += x;
					Cpx += x;
					break;	

				case P_ALBITE:
					ab += x;
					break;

				case P_ALMANDINE:
				case P_ANDRADITE:
				case P_GROSSULAR:
				case P_PYROPE:
				case P_SPESSARTINE:
					M += x;
					garnet += x;
					break;	
					
				case P_ANALCITE:
					F += x;
					anc += x;
					break;

				case P_ANDALUSITE:
				case P_KYANITE:
				case P_SILLIMANITE:
					M += x;
					anda += x;
					break;

				case P_ANNITE:	
				case P_NA_PHLOGOPITE:
				case P_PHLOGOPITE:
					M += x;
					bio += x;
					break;

				case P_ANORTHITE:
					an += x;
					P += x;
					break;

				case P_ANTHOPHYLLITE:
				case P_FE_ANTHOPHYLLITE:
					ho += x;
					M += x;
					break;

				case P_AKERMANITE:
				case P_BIXBYITE:
				case P_CLINOCHLORE:
				case P_DAPHNITE:
				case P_EPIDOTE:
				case P_Fe:
				case P_FE_STAUROLITE:
				case P_FERROCARPHOLITE:
				case P_GEHLENITE:
				case P_HEMATITE:
				case P_HEULANDITE:
				case P_HYDROXYAPATITE:
				case P_ILMENITE:
				case P_JADEITE:
				case P_MAGNESIOCARPHOLITE:
				case P_MAGNESIOFERRITE:
				case P_MAGNETITE:
				case P_P4:
				case P_PSEUDOWOLLASTONITE:
				case P_PYRRHOTITE:
				case P_PYROPHANITE:
				case P_PYROXMANGITE:
				case P_RHODONITE:
				case P_RUTILE:
				case P_SPHENE:
				case P_MG_STILPNOMELANE:
				case P_TALC:
				case P_TEPHROITE:
				case P_VESUVIANITE:
				case P_WOLLASTONITE:
				case P_WUSTITE:
					M += x;
					break;	

				case P_CLINOENSTATITE:
					hy += x;
					M += x;
					Cpx += x;
					break;		

				case P_CORDIERITE:
					cor += x;
					M += x;
					break;		

				case P_CORUNDUM:
					co += x;
					M += x;
					break;		

				case P_DIOPSIDE:
				case P_HEDENBERGITE:
					M += x;
					di += x;
					Cpx += x;
					break;

				case P_ENSTATITE:
				case P_PROTOENSTATITE:
					hy += x;
					M += x;
					Opx += x;
					break;		

				case P_FAYALITE:
				case P_FORSTERITE:
					M += x;
					ol += x;
					break;

				case P_HERCYNITE:
				case P_SPINEL:
				case P_ULVOSPINEL:
					M += x;
					spinel += x;
					break;	

				case P_KALSILITE:
				case P_LEUCITE:
					lc += x;
					F += x;
					break;

				case P_MICROCLINE:
				case P_SANIDINE:
					A += x;
					break;

				case P_MUSCOVITE:
					A += x;
					mc += x;
					break;

				case P_NEPHELINE:
					ne += x;
					F += x;
					break;

				case P_ALBITE_LIQUID:
				case P_ANORTHITE_LIQUID:
				case P_CORUNDUM_LIQUID:
				case P_DIOPSIDE_LIQUID:
				case P_O2:
				case P_SiO2_LIQUID:
				case P_H2:
				case P_H2O_LIQUID:
				case P_WATER_VAPOR:
				case P_Na2SiO3:
				case P_WOLLASTONITE_LIQUID:	
					break;

				case P_CRISTOBALITE:
				case P_QUARTZ:
				case P_TRIDYMITE:
					Q += x;
					break;

				default:
					cout << "need classify case for component " 
						<< state_phase[state_ph[i]].name << endl;
					exit(1);
			}
		}
	}

	An = an/(ab + an);
	if (An < 0.05)
		A += an + ab;
	else
		P += an + ab;

	double x1 = 0.0, x2 = 0.0, xm = 0.0;
	sort_mafic(di, "augite", x1, m1, x2, m2);
	sort_mafic(ol, "olivine", x1, m1, x2, m2);
	sort_mafic(anda, "andalusite", x1, m1, x2, m2);
	sort_mafic(bio, "biotite", x1, m1, x2, m2);
	sort_mafic(anc, "analcite", x1, m1, x2, m2);
	sort_mafic(cor, "cordierite", x1, m1, x2, m2);
	sort_mafic(hy, "hypersthene", x1, m1, x2, m2);
	sort_mafic(ho, "hornblende", x1, m1, x2, m2);
	sort_mafic(co, "corundum", x1, m1, x2, m2);
	sort_mafic(mc, "muscovite", x1, m1, x2, m2);

    sort_ultramafic(garnet, "garnet", xm, um);
	sort_ultramafic(spinel, "spinel", xm, um);

	if (xm<5) 
		um = um + "-bearing";
}

