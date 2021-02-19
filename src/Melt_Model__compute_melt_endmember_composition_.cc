/*
 * Melt_Model__compute_melt_endmember_composition_.cc
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

//-----------------------------------------------------------------------------//
/*! Calculate endmember composition of melt.
 * 
 * \param[in]  XM Molar composition of melt by chemical element
 * \param[out] x Endmember composition of melt
 */
void Melt_Model::compute_melt_endmember_composition_(double const XM[E_END],
                                                     double x[M_END]) const
{
	using namespace std;
	
	fill(x, x+M_END, 0.0);

	// Reorganize in analogy to CIPW norm into preferred melt phases.

	// 1: I have no phosphate melt

	// 2: I have no pyrite melt.

	// 3: I have no chromite melt.

	// 4: I hae no ilmenite melt.

	// 5: I have no fluorite melt.

	// 6: I have no calcite melt

	// 7: I have no zircon melt.

	// 8: Orthoclase

    double Q = 0.5*XM[E_K];
	x[M_KAlSi3O8] = Q;
	x[M_Al2O3] = 0.5*XM[E_AL] - Q;
	x[M_SiO2] = XM[E_SI] - 6*Q;

	// 9: Albite
	double Na2O = 0.5*XM[E_NA];
	Q = 0.5*min(Na2O, x[M_Al2O3]);
	x[M_NaAlSi3O8] = Q;
	Na2O -= 0.5*Q;
	x[M_Al2O3] -= 0.5*Q;
	x[M_SiO2] -= 3*Q;

	// 10: Anorthite
	x[M_CaO] = XM[E_CA];
	Q = min(x[M_CaO], x[M_Al2O3]);
	x[M_CaAl2Si2O8] = Q;
	x[M_CaO] -= Q;
	x[M_Al2O3] -= Q;
	x[M_SiO2] -= 2*Q;

	// 11: Acmite should come next, but I have no melt for it.

	// 12: Sodium metasilicate

	Q = Na2O;
	Na2O = 0.0;
	x[M_Na2SiO3] = Q;
	x[M_SiO2] -= Q;

	// 13: I have no melts for magnetite and am not presently implementing hematite.

	// 14-15: Magnesium diopside. I have no hedenbergite melt.

	x[M_MgO] = XM[E_MG];
	Q = min(x[M_MgO], x[M_CaO]);
	x[M_CaMgSi2O6] = Q;
	x[M_CaO] -= Q;
	x[M_MgO] -= Q;
	x[M_SiO2] -= 2*Q;

	// 16: Wollastonite

	Q = x[M_CaO];
	x[M_CaSiO3] = Q;
	x[M_CaO] = 0.0;
	x[M_SiO2] -= Q;

	// 17: Magnesium pyroxene (enstatite). I have no ferrosilite melt.

	Q = 0.5*x[M_MgO];
	x[M_Mg2Si2O6] = Q;
	x[M_MgO] = 0.0;
	x[M_SiO2] -= 2*Q;

	Q = 0.5*XM[E_FE];
	x[M_Fe2SiO4] = Q;
	x[M_SiO2] -= Q;

	if (x[M_SiO2]<0.0)
	{
		// Silica undersaturated

		double D = -x[M_SiO2];
		x[M_SiO2] = 0.0;

		// Enstatite to olivine.
		Q = min(D, x[M_Mg2Si2O6]);
		x[M_Mg2SiO4] = Q;
		x[M_Mg2Si2O6] -= Q;
		D -= Q;

		if (D>0.0)
		{
			// Albite to nepheline

			Q = min(0.5*D, x[M_NaAlSi3O8]);
			x[M_NaAlSiO4] += Q;
			x[M_NaAlSi3O8] -= Q;
			D -= 2*Q;

			if (D>0.0)
			{
				// Orthoclase to leucite

				Q = min(x[M_KAlSi3O8], D);
				x[M_KAlSi2O6] += Q;
				x[M_KAlSi3O8] -= Q;
				D -= Q;

				if (D>0.0)
				{
					// Do not have calcium orthosilicate. Go to lime instead.

					// Diopside to wollastonite and olivine.

					Q = min(x[M_CaMgSi2O6], 2*D);
					x[M_CaMgSi2O6] -= Q;
					x[M_CaSiO3] += Q;
					x[M_Mg2SiO4] += 0.5*Q;
					D -= 0.5*Q;


					if (D>0.0)
					{
						// Wollastonite to lime

						Q = min(D, x[M_CaSiO3]);
						x[M_CaO] += Q;
						x[M_CaSiO3] -= Q;
						D -= Q;

						if (D>0.0)
						{
							// Olivine to periclase

							Q = min(D, x[M_Mg2SiO4]);
							x[M_MgO] += 2*Q;
							x[M_Mg2SiO4] -= Q;
							D -= Q;

							// That's all I have melts for. If quartz is still deficient, this is not a valid melt.
							if (D<-1.0e-9)
							{
								Insist(false, "bad branch");
							}
						}
					}
				}
			}
		}
	}

	x[M_H2O] = 0.5*XM[E_H];

//	cout << "Melt CIPW:" << endl;
//	for (unsigned m=0; m<M_END; ++m)
//	{
//		if (x[m]>0.0)
//		{
//			cout << phase[endmember[m]].name << ": " << x[m] << endl;
//		}
//	}
}
