/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * iugs_classify.cc
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

#include "classify.hh"
#include "gui.hh"
#include "update.hh"

std::string mafic(std::string const &m2, std::string const &m1)
{
	if (m2.size()>0)
	{
		if (m1.size()>0)
		{
			return m2 + '-' + m1 + ' ';
		}
		else
		{
			return m2 + ' ';
		}
	}
	else
	{
		if (m1.size()>0)
		{
			return m1 + ' ';
		}
		else
		{
			return "";
		}
	}
}

//-----------------------------------------------------------------------------//
void IUGS_classify(State const &state)
{
	using namespace std;
	
	double Q, an, anc, ne, lc, An, A, P, F, M, Ol, Opx, Cpx;
	string m1, m2, um;

	classify_components(state, Q, an, anc, ne, lc, An, A, P, F, M, Ol, Opx, Cpx,
	                    m1, m2, um);

	double Pfrac = 100*P/(A+P);

	if (M<90)
	{
		if (F==0)
		{
			if (Q<=5)
			{
				if (Pfrac<=10)
				{
					double crit = 23 + Pfrac*17./65.;
					if (M>crit)
						text_iugs->set_text(mafic(m2, m1) + "alkali-feldspar melasyenite");
					else
						text_iugs->set_text("Alkali-feldspar syenite");
				}
				else if (Pfrac<=35)
				{
					text_iugs->set_text("SYENITE\ntrachyte");
				}
				else if (Pfrac<=65)
				{
					text_iugs->set_text(mafic(m2, m1) + "monzonite");
				}
				else if (Pfrac<=90)
				{
					if (An>50)
					{
						text_iugs->set_text("MONZO-GABBRO\nbasalt");
					}
					else
					{
						text_iugs->set_text(mafic(m2, m1) + "monzodiorite");
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
						text_iugs->set_text(mafic(m2, m1) + "diorite");
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
						text_iugs->set_text(mafic(m2,m1)+"quartz monzodiorite");
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
						text_iugs->set_text(mafic(m2, m1) + "quartz diorite");
					}
				}
			}
			else if (Q<=60)
			{
				if (Pfrac<=10)
				{
					if (M>20)				
						text_iugs->set_text("Alkali-feldspar melagranite");
					else
						text_iugs->set_text(mafic(m2, m1) + "alkali-feldspar granite");
				}
				else if (Pfrac<=65)
				{
					text_iugs->set_text(mafic(m2, m1) + "granite");
				}
				else if (Pfrac<=90)
				{
					text_iugs->set_text(mafic(m2, m1) + "granodiorite");
				}
				else
				{
					text_iugs->set_text(mafic(m2, m1) + "tonalite");	
				}			
			}
			else if (Q<=90)
			{
				text_iugs->set_text("QUARTZ-RICH GRANITOID\nrhyolite");			
			}
			else
			{
				text_iugs->set_text(mafic(m2, m1) + "quartzolite");
			}
		}
		else // F>0
		{
			if (F<=10)
			{
				if (Pfrac<=10)
				{
					text_iugs->set_text(mafic(m2, m1) + "nepheline-bearing alkali-feldspar syenite");
				}
				else if (Pfrac<=35)
				{
					text_iugs->set_text("Nepheline-bearing syenite");	
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
						text_iugs->set_text(mafic(m2, m1) + " nepheline-bearing monzodiorite");	
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
						if (lc>ne)
						text_iugs->set_text(mafic(m2, m1) + "leucite-bearing diorite");	
						else
						text_iugs->set_text(mafic(m2, m1) + "nepheline-bearing diorite");	
					}
				}
			}
			else if (F<=60)
			{
				if (Pfrac<=10)
				{
					if (M>40)
					{
						if (ne>lc)
							text_iugs->set_text(mafic(m2, m1) + "nepheline melasyenite");	
						else
							text_iugs->set_text(mafic(m2, m1) + "leucite melasyenite");	
					}
					else
					{
						if (ne>lc)
							text_iugs->set_text(mafic(m2, m1) + "nepheline syenite");	
						else
							text_iugs->set_text(mafic(m2, m1) + "leucite syenite");	
					}
				}
				else if (Pfrac<=50)
				{
					if (ne>lc)
				  	  text_iugs->set_text(mafic(m2, m1) + "nepheline monzosyenite");
					else
				  	  text_iugs->set_text(mafic(m2, m1) + "leucite monzosyenite");
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
						if (ne>lc)
						  text_iugs->set_text(mafic(m2, m1) + "nepheline diorite");	
						else
						  text_iugs->set_text(mafic(m2, m1) + "leucite diorite");	
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
		double norm = Ol + Opx + Cpx;
		if (norm>0)
		{
			norm = 100/norm;
			Ol *= norm;
			Opx *= norm;
			Cpx *= norm;
			if (Ol>90)
			{
				text_iugs->set_text("<ultramafic>");	
			}
			else if (Ol > 40)
			{
				if (Cpx<5)
				{
					text_iugs->set_text("harzburgite");	
				}
				else if (Opx<5)
				{
					text_iugs->set_text(um + " wehrlite");	
				}
				else
				{
					text_iugs->set_text(um + " lherzolite");	
				}
			}
			else 
			{
				text_iugs->set_text("<ultramafic>");	
			}
		}
		else
		{
			text_iugs->set_text("<ultramafic>");	
		}
	}
}

