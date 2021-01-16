/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * update.cc
 * Copyright (C) 2021 Kent G. Budge <kgb@kgbudge.com>
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
void update()
{
	using namespace std;
	
	read_weight_composition();

	calculate_molar_composition();
	if (ntotal==0.0) return;

	populate_mole_composition_box();

	// Magnesium number

	text_MgNo->set_text(tostring((int)(1000*nMgO/(nMgO+nFeO)+0.5)/10.));

	try
	{
		calculate_CIPW_norm();
		populate_CIPW_norm_box();
		populate_cation_norm_box();
		
		text_DI->set_text(tostring(calculate_differentiation_index(), 2));
		
		text_tas->set_text(TAS_classify());

		text_ibc_family->set_text(IBC_classify());

		// Determine IUGS classification

		// Convert to volume fraction

		ap /=3.16;
		il /=4.9;
		orc /=2.63;
		ab /=2.6;
		an /=2.75;
		nAl2O3 /=3.95;
		ac /=3.55;
		nNa2O /=2.61;
		mt /=5.18;
		nFe2O3 /=5.26;
		di /=(3.28*wMgFe + 3.56*(1-wMgFe));
		nCaO /=2.86;
		hy /=(3.4*wMgFe + 3.9*(1-wMgFe));
		nSiO2 /=2.65;
		ol /=(3.2*wMgFe + 4.5*(1-wMgFe));
		tn /=3.48;
		nTiO2 /=4.23;
		pf /=3.98;
		ne /=2.59;
		nK2O /=2.61; // lacking density of potassium silicate, we assume same molar volume as sodium metasilicate
		lc /=2.50;
		cs /=2.9;
		pr /=5.0;
		cm /=4.8;
		kp /=2.59;
		ph /=2.39;
		kk /=2.32; 

		double total = ap + il + orc + ab + an + nAl2O3 + ac + nNa2O + mt + nFe2O3 + di
			+ nCaO + hy + nSiO2 + ol + tn + nTiO2 + pf
			+ ne + nK2O + lc + cs + pr + cm + kp + ph + kk;

		double rnorm = 100.0/total;

		ap *= rnorm;
		il *= rnorm;
		orc *= rnorm;
		ab *= rnorm;
		an *= rnorm;
		nAl2O3 *= rnorm;
		ac *= rnorm;
		nNa2O *= rnorm;
		mt *= rnorm;
		nFe2O3 *= rnorm;
		di *= rnorm;
		nCaO *= rnorm;
		hy *= rnorm;
		nSiO2 *= rnorm;
		ol *= rnorm;
		tn *= rnorm;
		nTiO2 *= rnorm;
		pf *= rnorm;
		ne *= rnorm;
		nK2O *= rnorm;
		lc *= rnorm;
		cs *= rnorm;
		pr *= rnorm;
		cm *= rnorm;
		kp *= rnorm;
		ph *= rnorm;
		kk *= rnorm;

		double Q = nSiO2;
		double An = 100*an/(an+ab);
		double A = orc + (An>0.05? 0 : ab);
		double P = an  + (An>0.05? ab : 0);
		double F = ne + lc + kp;

		double M = ol + ac + di + hy + mt + cm + nFe2O3 + il + tn + pf + nTiO2
			+ pr + ap + nAl2O3;

		double Pfrac = 100*P/(A+P);

		if (M<90)
		{
			if (F==0)
			{
				if (Q<=5)
				{
					if (Pfrac<=10)
					{
						text_iugs->set_text("ALKALI-FELDSPAR SYENITE\ntrachyte");
					}
					else if (Pfrac<=35)
					{
						text_iugs->set_text("SYENITE\ntrachyte");
					}
					else if (Pfrac<=65)
					{
						text_iugs->set_text("MONZONITE\nlatite");
					}
					else if (Pfrac<=90)
					{
						if (An>50)
						{
							text_iugs->set_text("MONZO-GABBRO\nbasalt");
						}
						else
						{
							text_iugs->set_text("MONZO-DIORITE\nandesite");
						}
					}
					else
					{
						if (An>50)
						{
							text_iugs->set_text("GABBRO\nbasalt");
						}
						else
						{
							text_iugs->set_text("DIORITE\nandesite");
						}
					}
				}
				else if (Q<=20)
				{
					if (Pfrac<=10)
					{
						text_iugs->set_text("QUARTZ ALKALI-FELDSPAR SYENITE\ntrachyte");
					}
					else if (Pfrac<=35)
					{
						text_iugs->set_text("QUARTZ SYENITE\nquartz trachyte");
					}
					else if (Pfrac<=65)
					{
						text_iugs->set_text("QUARTZ MONZONITE\nquartz latite");
					}
					else if (Pfrac<=90)
					{
						if (An>50)
						{
							text_iugs->set_text("QUARTZ MONZOGABBRO\nandesite");
						}
						else
						{
							text_iugs->set_text("QUARTZ MONZODIORITE\nandesite");
						}
					}
					else
					{
						if (an>90)
						{
							text_iugs->set_text("QUARTZ ANORTHOSITE");
						}
						else if (An>50)
						{
							text_iugs->set_text("QUARTZ GABBRO\nbasalt");
						}
						else
						{
							text_iugs->set_text("QUARTZ DIORITE\nandesite");
						}
					}
				}
				else if (Q<=60)
				{
					if (Pfrac<=10)
					{
						text_iugs->set_text("ALKALI-FELDSPAR GRANITE\nalkali-feldspar rhyolite");
					}
					else if (Pfrac<=65)
					{
						text_iugs->set_text("GRANITE\nrhyolite");
					}
					else if (Pfrac<=90)
					{
						text_iugs->set_text("GRANODIORITE\ndacite");
					}
					else
					{
						text_iugs->set_text("TONALITE\ndacite");	
					}			
				}
				else if (Q<=90)
				{
					text_iugs->set_text("QUARTZ-RICH GRANITOID\nrhyolite");			
				}
				else
				{
					text_iugs->set_text("QUARTZOLITE\nrhyolite");
				}
			}
			else // F>0
			{
				if (F<=10)
				{
					if (Pfrac<=10)
					{
						text_iugs->set_text("FELDSPATHOID-BEARING ALKALI-FELDSPAR SYENITE\nfeldspathoid-bearing trachyte");
					}
					else if (Pfrac<=35)
					{
						text_iugs->set_text("FELDSPATHOID-BEARING SYENITE\nfeldspathoid-bearing trachyte");	
					}
					else if (Pfrac<=65)
					{
						text_iugs->set_text("FELDSPATHOID-BEARING MONZONITE\nfeldspathoid-bearing latite");	
					}
					else if (Pfrac<=90)
					{
						if (An>50)
						{
							text_iugs->set_text("FELDSPATHOID-BEARING MONZOGABBRO\nfeldspathoid-bearing basalt");	
						}
						else
						{
							text_iugs->set_text("FELDSPATHOID-BEARING MONZODIORITE\nfeldspathoid-bearing basalt");	
						}
					}
					else
					{
						if (An>50)
						{
							text_iugs->set_text("FELDSPATHOID-BEARING GABBRO\nfeldspathoid-bearing basalt");	
						}
						else
						{
							text_iugs->set_text("FELDSPATHOID-BEARING DIORITE\nfeldspathoid-bearing basalt");	
						}
					}
				}
				else if (F<=60)
				{
					if (Pfrac<=10)
					{
						text_iugs->set_text("FELDSPATHOID SYENITE\nphonolite");	
					}
					else if (Pfrac<=50)
					{
						text_iugs->set_text("FELDSPATHOID MONZOSYENITE\ntephritic phonolite");	
					}
					else if (Pfrac<=90)
					{
						if (An>50)
						{
							text_iugs->set_text("FELDSPATHOID MONZOGABBRO\nphonolitic tephrite");	
						}
						else
						{
							text_iugs->set_text("FELDSPATHOID MONZODIORITE\nphonolitic tephrite");	
						}
					}
					else
					{
						if (An>50)
						{
							text_iugs->set_text("FELDSPATHOID GABBRO\nbasanite");	
						}
						else
						{
							text_iugs->set_text("FELDSPATHOID DIORITE\ntephrite");	
						}
					}
				}
				else if (F<=90)
				{
					if (Pfrac<=50)
					{
						text_iugs->set_text("PHONOLITIC FOIDITE\nfoidite");						
					}
					else
					{
						text_iugs->set_text("TEPHRITIC FOIDITE\nfoidite");						
					}
				}
				else
				{
					text_iugs->set_text("FOIDITE\nfoidite");	
				}
			}
		}
		else
		{
			text_iugs->set_text("<ultramafic>");			
		}
	}
	catch (bad_CIPW const &err)
	{
		Gtk::MessageDialog dialog(*main_win, "CIPW Norm Fail");
		dialog.set_secondary_text(err.what());
		dialog.run();

		clear_CIPW_norm_box();
		clear_cation_norm_box();
	}
}
