/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * gui.hh
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

#ifndef melt_gui_hh
#define melt_gui_hh

#include <gtkmm.h>

#ifndef EXTERN
#define EXTERN extern
#endif

//-----------------------------------------------------------------------------//
class CIPW_Columns : public Gtk::TreeModelColumnRecord
{
public:

  CIPW_Columns()
    { add(m_col_text); add(m_col_number); }

  Gtk::TreeModelColumn<Glib::ustring> m_col_text;
  Gtk::TreeModelColumn<double> m_col_number;
};

//-----------------------------------------------------------------------------//

EXTERN Gtk::Window* main_win;

EXTERN Gtk::Entry *entry_name;

EXTERN Gtk::Entry *entry_SiO2;
EXTERN Gtk::Entry *entry_TiO2;
EXTERN Gtk::Entry *entry_Al2O3;
EXTERN Gtk::Entry *entry_Fe2O3;
EXTERN Gtk::Entry *entry_FeO;
EXTERN Gtk::Entry *entry_MnO;
EXTERN Gtk::Entry *entry_MgO;
EXTERN Gtk::Entry *entry_CaO;
EXTERN Gtk::Entry *entry_Na2O;
EXTERN Gtk::Entry *entry_K2O;
EXTERN Gtk::Entry *entry_P2O5;
EXTERN Gtk::Entry *entry_S;
EXTERN Gtk::Entry *entry_Cr2O3;
EXTERN Gtk::Entry *entry_ZrO2;
EXTERN Gtk::Entry *entry_H2O;
EXTERN Gtk::Entry *entry_CO2;
EXTERN Gtk::Entry *entry_Cl;
EXTERN Gtk::Entry *entry_T;
EXTERN Gtk::Entry *entry_P;
EXTERN Gtk::Entry *entry_over_rho;
EXTERN Gtk::Entry *entry_depth;


EXTERN Gtk::Label *text_nSiO2;
EXTERN Gtk::Label *text_nTiO2;
EXTERN Gtk::Label *text_nAl2O3;
EXTERN Gtk::Label *text_nFe2O3;
EXTERN Gtk::Label *text_nFeO;
EXTERN Gtk::Label *text_nMnO;
EXTERN Gtk::Label *text_nMgO;
EXTERN Gtk::Label *text_nCaO;
EXTERN Gtk::Label *text_nNa2O;
EXTERN Gtk::Label *text_nK2O;
EXTERN Gtk::Label *text_nP2O5;
EXTERN Gtk::Label *text_nS;
EXTERN Gtk::Label *text_nCr2O3;
EXTERN Gtk::Label *text_nZrO2;
EXTERN Gtk::Label *text_nH2O;
EXTERN Gtk::Label *text_nCO2;
EXTERN Gtk::Label *text_nCl;

EXTERN Gtk::Label *text_total;

EXTERN Gtk::Label *text_iugs;

EXTERN Gtk::RadioButton *button_byweight;
EXTERN Gtk::RadioButton *button_oxygen_by_composition;
EXTERN Gtk::RadioButton *button_oxygen_specified;
EXTERN Gtk::RadioButton *button_oxygen_FMQ;
EXTERN Gtk::RadioButton *button_molar_melt;
EXTERN Gtk::Entry *entry_pO2;

EXTERN Gtk::TreeView *tree_view_CIPW;

EXTERN std::string file;
EXTERN std::string name;

//-----------------------------------------------------------------------------//

void on_activate();
void on_activate_depth();
void on_activate_P();
void on_new();
void on_normalize();
void on_open();
void on_quit();
void on_save();
void on_save_as();

void do_save();
void update();

#endif
