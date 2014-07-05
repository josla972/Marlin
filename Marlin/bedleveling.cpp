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
#include "Marlin.h"
#include "bedleveling.h"
#include "planner.h"
#include "stepper.h"
#include "language.h"

#define DEFINE_PGM_READ_ANY(type, reader)       \
    static inline type pgm_read_any(const type *p)  \
    { return pgm_read_##reader##_near(p); }

DEFINE_PGM_READ_ANY(float,       float);
DEFINE_PGM_READ_ANY(signed char, byte);

#define XYZ_CONSTS_FROM_CONFIG(type, array, CONFIG) \
static const PROGMEM type array##_P[3] =        \
    { X_##CONFIG, Y_##CONFIG, Z_##CONFIG };     \
static inline type array(int axis)          \
    { return pgm_read_any(&array##_P[axis]); }

static float feedrate;
static int probePointCounter;
XYZ_CONSTS_FROM_CONFIG(float, home_retract_mm, HOME_RETRACT_MM);

#define SQUARE(x) (x*x) // Parantheses used to keep priority in case of A/SQUARE(B) = A/(B*B) and not A/B*B
void probeLine(int start_x, int start_y, int xStep, int yStep, double *eqnAMatrix, double *eqnBVector);
void probeEdges(double *eqnAMatrix, double *eqnBVector);
void engage_z_probe();
void retract_z_probe();
void run_z_probe();

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
// / Josef Larsson - josla972@student.liu.se
double calculate_z_variance(int points, double *plane_equation_coefficents, double *eqnAMatrix, double *eqnBVector)
{
	double MSE = 0; // Mean squared error
	double bias = 0;
	for (int i = 0; i < points; i++)
	{ 
		double error = eqnBVector[i] - (plane_equation_coefficents[0]*eqnAMatrix[i + 0*points] + plane_equation_coefficents[1]*eqnAMatrix[i + 1*points] + plane_equation_coefficents[2]); // e = z_measured - z = z_measured - (ax + by + d)
		MSE += SQUARE(error)/(float)points;
		bias += error/(float)points;
  		SERIAL_PROTOCOLPGM(" x: ");
  		SERIAL_PROTOCOL(eqnAMatrix[i + 0*points]);
  		SERIAL_PROTOCOLPGM(" y: ");
  		SERIAL_PROTOCOL(eqnAMatrix[i + 1*points]);
  		SERIAL_PROTOCOLPGM(" z: ");
  		SERIAL_PROTOCOL(eqnBVector[i]);
  		SERIAL_PROTOCOLPGM("\n");
	}
	if (MSE - SQUARE(bias) < 0) 
	{
		SERIAL_PROTOCOLPGM("Shit hit the fan: Negative variance!");
	}
	return MSE - SQUARE(bias); 
}

void probePoint(int xProbe, int yProbe, double *eqnAMatrix, double *eqnBVector)
{
	float z_before = (probePointCounter == 0 ? Z_RAISE_BEFORE_PROBING : current_position[Z_AXIS] + Z_RAISE_BETWEEN_PROBINGS);
	
	float measured_z = probe_pt(xProbe, yProbe, z_before);
	
	eqnBVector[probePointCounter] = measured_z;

#ifdef AUTO_BED_LEVELING_EDGES
 	SERIAL_PROTOCOLPGM(MSG_BED);
 	SERIAL_PROTOCOLPGM("Filling eqnAMatrix[");
 	SERIAL_PROTOCOL(probePointCounter + 0*(AUTO_BED_LEVELING_GRID_POINTS-1)*4);
 	SERIAL_PROTOCOLPGM("], with: ");
 	SERIAL_PROTOCOL(xProbe);
 	SERIAL_PROTOCOLPGM("\n and");
 	SERIAL_PROTOCOLPGM("Filling eqnAMatrix[");
 	SERIAL_PROTOCOL(probePointCounter + 1*(AUTO_BED_LEVELING_GRID_POINTS-1)*4);
 	SERIAL_PROTOCOLPGM("], with: ");
 	SERIAL_PROTOCOL(yProbe);
 	SERIAL_PROTOCOLPGM("\n");
	eqnAMatrix[probePointCounter + 0*(AUTO_BED_LEVELING_GRID_POINTS-1)*4] = xProbe;
	eqnAMatrix[probePointCounter + 1*(AUTO_BED_LEVELING_GRID_POINTS-1)*4] = yProbe;
	eqnAMatrix[probePointCounter + 2*(AUTO_BED_LEVELING_GRID_POINTS-1)*4] = 1;
#else
	eqnAMatrix[probePointCounter + 0*AUTO_BED_LEVELING_GRID_POINTS*AUTO_BED_LEVELING_GRID_POINTS] = xProbe;
	eqnAMatrix[probePointCounter + 1*AUTO_BED_LEVELING_GRID_POINTS*AUTO_BED_LEVELING_GRID_POINTS] = yProbe;
	eqnAMatrix[probePointCounter + 2*AUTO_BED_LEVELING_GRID_POINTS*AUTO_BED_LEVELING_GRID_POINTS] = 1;
#endif
}

void probeEdges(double *eqnAMatrix, double *eqnBVector) 
{
	int xGridSpacing = (X_MAX_POS - X_MIN_POS) / (AUTO_BED_LEVELING_GRID_POINTS-1);
	int yGridSpacing = (Y_MAX_POS - Y_MIN_POS) / (AUTO_BED_LEVELING_GRID_POINTS-1);
	probePointCounter = 0;
	probeLine(X_MIN_POS, Y_MIN_POS, xGridSpacing, 0, eqnAMatrix, eqnBVector); 
	probeLine(X_MAX_POS, Y_MIN_POS, 0, yGridSpacing, eqnAMatrix, eqnBVector); 
	probeLine(X_MAX_POS, Y_MAX_POS, -xGridSpacing, 0, eqnAMatrix, eqnBVector); 
	probeLine(X_MIN_POS, Y_MAX_POS, 0, -yGridSpacing, eqnAMatrix, eqnBVector); 
}

void probeLine(int start_x, int start_y, int xStep, int yStep, double *eqnAMatrix, double *eqnBVector)
{
	for (int i=0, xProbe = start_x, yProbe = start_y; i < AUTO_BED_LEVELING_GRID_POINTS-1; i++)
	{
		probePoint(xProbe, yProbe, eqnAMatrix, eqnBVector);	
		xProbe += xStep;
		yProbe += yStep;
		probePointCounter++;
	}
}

void do_blocking_move_to(float x, float y, float z) {
//    float oldFeedRate = feedrate;

 //   feedrate = new_feedrate;

    current_position[X_AXIS] = x;
    current_position[Y_AXIS] = y;
    current_position[Z_AXIS] = z;
    plan_buffer_line(current_position[X_AXIS], current_position[Y_AXIS], current_position[Z_AXIS], current_position[E_AXIS], feedrate/60, active_extruder);
    st_synchronize();

//    feedrate = oldFeedRate;
}

void do_blocking_move_relative(float offset_x, float offset_y, float offset_z) {
    do_blocking_move_to(current_position[X_AXIS] + offset_x, current_position[Y_AXIS] + offset_y, current_position[Z_AXIS] + offset_z, feedrate);
}

void set_feedrate(float new_feedrate)
{
	feedrate = new_feedrate;
}

// Probe bed height at position (x,y), returns the measured z value
float probe_pt(float x, float y, float z_before) {
  float old_feedrate = feedrate;
  // move to right place
  do_blocking_move_to(current_position[X_AXIS], current_position[Y_AXIS], z_before);
  feedrate = XY_TRAVEL_SPEED;
  do_blocking_move_to(x - X_PROBE_OFFSET_FROM_EXTRUDER, y - Y_PROBE_OFFSET_FROM_EXTRUDER, current_position[Z_AXIS]);
  feedrate = old_feedrate;

  engage_z_probe();   // Engage Z Servo endstop if available
  run_z_probe();
  float measured_z = current_position[Z_AXIS];
  retract_z_probe();

  SERIAL_PROTOCOLPGM(MSG_BED);
  SERIAL_PROTOCOLPGM(" x: ");
  SERIAL_PROTOCOL(x);
  SERIAL_PROTOCOLPGM(" y: ");
  SERIAL_PROTOCOL(y);
  SERIAL_PROTOCOLPGM(" z: ");
  SERIAL_PROTOCOL(measured_z);
  SERIAL_PROTOCOLPGM("\n");
  return measured_z;
}

void run_z_probe() 
{
    plan_bed_level_matrix.set_to_identity();
    feedrate = homing_feedrate[Z_AXIS];

    // move down until you find the bed
    float zPosition = -10;
    plan_buffer_line(current_position[X_AXIS], current_position[Y_AXIS], zPosition, current_position[E_AXIS], feedrate/60, active_extruder);
    st_synchronize();

    // we have to let the planner know where we are right now as it is not where we said to go.
    zPosition = st_get_position_mm(Z_AXIS);
    plan_set_position(current_position[X_AXIS], current_position[Y_AXIS], zPosition, current_position[E_AXIS]);

    // move up the retract distance
    zPosition += home_retract_mm(Z_AXIS);
    plan_buffer_line(current_position[X_AXIS], current_position[Y_AXIS], zPosition, current_position[E_AXIS], feedrate/60, active_extruder);
    st_synchronize();

    // move back down slowly to find bed
    feedrate = homing_feedrate[Z_AXIS]/4;
    zPosition -= home_retract_mm(Z_AXIS) * 2;
    plan_buffer_line(current_position[X_AXIS], current_position[Y_AXIS], zPosition, current_position[E_AXIS], feedrate/60, active_extruder);
    st_synchronize();

    current_position[Z_AXIS] = st_get_position_mm(Z_AXIS);
    // make sure the planner knows where we are as it may be a bit different than we last said to move to
    plan_set_position(current_position[X_AXIS], current_position[Y_AXIS], current_position[Z_AXIS], current_position[E_AXIS]);
}

