/*
 * Melt_Model.cc
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

#include "Melt_Model.hh"

#include "Model.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

//-----------------------------------------------------------------------------//
/*! Create a Melt_Model reflecting possible melting of an existing State.
 *
 * \param state Current state of the mineral ensemble.
 */ 
 Melt_Model::Melt_Model(State const &state)
:
 T_(state.T()), P_(state.P()), phase_(state.phase()), Gf_(state.Gf()), 
 is_fusible_(state.Gf().size())
{
	 using namespace std;
	 
	// Compute the fully melted composition (excluding phases for which we
	// do not have a liquid counterpart). Determine free energy contribution
	// of these nonfusible phases.
	 
	Gfr_ = 0.0;  // Gibbs free energy contribution of nonfusible phases 
    fill(Z_, Z_+E_END, 0.0); // Prepare to accumulate sample composition.
	                         //This will also be the starting melt. 
	 
	cout << "Sample phases:" << endl;
	double const *const state_X = state.X();
	int const *const state_ph = state.ph();
	unsigned const NPL = phase_.size();
	for (unsigned i=0; i<E_END; ++i)
	{
		double const x = state_X[i];
		double xpr[E_END];  // to store composition of phase
		fill(xpr, xpr+E_END, 0.0); 
		if (x>0.0)  // Is this phase actually present?
		{
			unsigned p = state_ph[i];
			Phase const &ph = phase_[p];	
			cout << ph.name << endl;
			unsigned const N = ph.nz; // Number of elements in the phase
			double xO = 0.0;          // To accumulate oxygen balance of the phase
			bool fusible = true;
			for (unsigned j=0; j<N && fusible; ++j) 
			{
				double const xj = x*ph.n[j];  // Number of moles of element j in the phase.
				switch(ph.z[j])             // Switch on the element atomic number
				{
					case E_H:
						xO -= 0.5*xj;
						xpr[E_H] += xj;
						break; 

				    case E_C:
						xO -= 2*xj;
						xpr[E_C] += xj;
						break;
						
					case E_O:
						xpr[E_O] += xj;
						xO += xj;
						break;

					case E_NA:
						xO -= 0.5*xj;
						xpr[E_NA] += xj;
						break;

					case E_MG:
						xO -= xj;
						xpr[E_MG] += xj;
						break;

					case E_AL:
						xpr[E_AL] += xj;
						xO -= 1.5*xj;
						break;

					case E_SI:
						xpr[E_SI] += xj;
						xO -= 2*xj;
						break;

					case E_S:
						xpr[E_S] += xj;
						break;

					case E_CL:
						xpr[E_CL] += xj;
						xO += 0.5*xj;
						break;

					case E_K:
						xpr[E_K] += xj;
						xO -= 0.5*xj;
						break;

					case E_CA:
						xO -= xj;
						xpr[E_CA] += xj;
						break;

					case E_FE:
						xpr[E_FE] += xj;
						break;

					default:
						// non-fusible mineral
						cout << "phase " << ph.name << " cannot melt." << endl;
						fusible = false;
						break;
				}
			}
			if (!fusible || fabs(xO-xpr[E_FE])>1e-9) // at present, cannot handle ferric or oxidized sulfur melts
			{
				cout << "  Not fusible" << endl;
				Gfr_ += x*Gf_[p];
			}
			else
			{
				for (unsigned m=0; m<E_END; ++m)
				{
					Z_[m] += xpr[m];
				}
			}
		}
	}
	cout << "Full melt elemental molar composition:" << endl << defaultfloat;
	cnorm_ = 0.0;
	for (unsigned i=0; i<E_END; ++i)
	{
		Check(Z_[i]>=0.0);
		cnorm_ += Z_[i];
		if (Z_[i]>0)
		{
	  	  cout << element_name[i] << ": " << setprecision(3) << Z_[i] << endl;
		}
	}

	// Find all potentiall fusible phases
	 bool const *const state_is_element_active = state.is_element_active();
	 for (unsigned p=0; p<NPL; ++p)
	 {
		 Phase const &ph = phase_[p];	
		 cout << ph.name << endl;
		 unsigned const N = ph.nz; // Number of elements in the phase
		 double xO = 0.0;          // To accumulate oxygen balance of the phase
		 double xFe = 0.0;
		 bool fusible = true;
		 for (unsigned j=0; j<N && fusible; ++j) 
		 {
			 double const xj = ph.n[j];  // Number of moles of element j in the phase.
			 unsigned const z = ph.z[j];
			 fusible = fusible && (state_is_element_active[z]);
			 switch(z)             // Switch on the element atomic number
			 {
				 case E_H:
					 xO -= 0.5*xj;
					 break; 

				 case E_C:
					 xO -= 2*xj;
					 break;

				 case E_O:
					 xO += xj;
					 break;

				 case E_NA:
					 xO -= 0.5*xj;
					 break;

				 case E_MG:
					 xO -= xj;
					 break;

				 case E_AL:
					 xO -= 1.5*xj;
					 break;

				 case E_SI:
					 xO -= 2*xj;
					 break;

				 case E_S:
					 break;

				 case E_CL:
					 xO += 0.5*xj;
					 break;

				 case E_K:
					 xO -= 0.5*xj;
					 break;

				 case E_CA:
					 xO -= xj;
					 break;

				 case E_FE:
					 xFe += xj;
					 break;

				 default:
					 // non-fusible mineral
					 cout << "phase " << ph.name << " cannot melt." << endl;
					 fusible = false;
					 break;
			 }
			 is_fusible_[p] = fusible && fabs(xO-xFe)<1e-9;
			 // at present, cannot handle ferric or oxidized sulfur melts
		 }
	 }

	 // Compute end-member melt phase free energies
	 for (unsigned i=0; i<M_END; ++i)
	 {
		 Phase const &ph = ::phase[melt_endmember[i]];
		 Gfm_[i] = ph.model->Gf(ph, T_, P_);
	 }
 }

//-----------------------------------------------------------------------------//
double Melt_Model::Gf(double const XP[P_END]) const
{
	using namespace std;
	
	unsigned const NM = NP_;

	double xm[E_END];
	for (unsigned i=0; i<E_END; ++i)
	{
		xm[i] = Z_[i];
	}
	// Modify for amount of each fusible phase crystallized out
	double Result = Gfr_;
	for (unsigned i=0; i<NM; ++i)
	{
		double const x = XP[i];
		if (x>0.0)
		{
			Phase const &phase = phase_[i];  
			unsigned const N = phase.nz;
			for (unsigned j=0; j<N; ++j)
			{
				unsigned z = phase.z[j];
				xm[z] -= x*phase.n[j];
				xm[z] = max(0.0, xm[z]);
			}
			Result += x*Gf_[i];
		}
	}
    Result += this->Gfm(xm);

	return Result;
}

char const * const endmember_element_name[] =
{
		"H",
		"Si",
		"Al",
		"Mg",
		"Fe(+2)",
		"Ca",
		"Na",
		"K"
};

