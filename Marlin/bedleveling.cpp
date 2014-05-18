/* This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */
#include "Configuration.h"

#define SQUARE(x) (x*x) // Parantheses used to keep priority in case of A/SQUARE(B) = A/(B*B) and not A/B*B

// I think we should move parts or all bedleveling code here. Marlin_main.cpp is getting messy!

// Calculate the variance of the error of our LS estimate using the equation:
// z = ax + by + d
//
// let the error e = z_measured - z 
// var e = E[(e- E[e])^2] = ... = / Wikipedia, variance/ = E[e^2] - (E[e])^2 = MSE - bias^2 
// Where MSE is an abbreviation for Mean squared error. 
//
// If the square root of our variance is bigger than acceptable, we might want to rerun the 
// estimation algorithm to ensure that we do not dig in to the print surface, or go too high so that nothing adheres to the print surface.
// An estimate is useless if we do not know how good it is!
double calculate_z_variance(double *plane_equation_coefficents, double *eqnAMatrix, double *eqnBVector)
{
	double MSE = 0; // Mean squared error
	double bias = 0;
	for (int i = 0; i < AUTO_BED_LEVELING_GRID_POINTS; i++)
	{ 
		double error = eqnBVector[i] - (plane_equation_coefficents[0]*eqnAMatrix[i + 0*SQUARE(AUTO_BED_LEVELING_GRID_POINTS)] + plane_equation_coefficents[1]*eqnAMatrix[i + 1*SQUARE(AUTO_BED_LEVELING_GRID_POINTS)] + plane_equation_coefficents[2]); // e = z_measured - z = z_measured - (ax + by + d)
		MSE += SQUARE(error)/AUTO_BED_LEVELING_GRID_POINTS;
		bias += error/AUTO_BED_LEVELING_GRID_POINTS;
	}
	if (MSE - SQUARE(bias) < 0) 
	{
		SERIAL_PROTOCOLPGM("Shit hit the fan: Negative variance!");
	}
	return MSE - SQUARE(bias); 
}
