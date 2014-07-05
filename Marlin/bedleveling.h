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
double calculate_z_variance(int points, double *plane_equation_coefficents, double *eqnAMatrix, double *eqnBVector);
void do_blocking_move_to(float x, float y, float z, float feedrate);
void do_blocking_move_relative(float offset_x, float offset_y, float offset_z);
float probe_pt(float x, float y, float z_before);
void set_feedrate(float new_feedrate);
void probeEdges(double *eqnAMatrix, double *eqnBVector);
// TODO: Perhaps we want to make these static again.
void engage_z_probe();
void retract_z_probe();
void run_z_probe();
