/*
 * melt.cc
 * Copyright (C) 2019 Kent G. Budge <kgb@kgbudge.com>
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

/*
#include <numeric>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_sort_vector.h"

#include "ds++/Assert.hh"

#include "algorithm.hh"
#include "constants.hh"
#include "D1.hh"
#include "element.hh"
#include "Melt_Model.hh"
#include "Model.hh"
#include "Phase.hh"
#include "State.hh"
 */

#include "Melt_Model.hh"

//-----------------------------------------------------------------------------//
void Melt_Model::compute_current_melt_composition_(double const XP[], double xm[]) const
{
	for (unsigned i=0; i<E_END; ++i)
	{
		xm[i] = Z_[i];
	}
	for (unsigned i=0; i<NP_; ++i)
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
			}
		}
	}
}
