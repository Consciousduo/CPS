/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_
#define DISTANCE_FACTOR 5.0


void computeAcceleration(struct particleSystem * cps, struct point a[N_particle]);

void printMatrixC(int row, int column, double matrix[][14]);
void printMatrixT(int row, int column, double matrix[][9]);
void printMatrixA(int row, int column, double matrix[][23]);




// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct particleSystem * cps);
#endif

