/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_
#define DISTANCE_FACTOR 5.0

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

void computeAcceleration(struct particleSystem * cps, struct point a[N_particle]);


// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct particleSystem * cps);
#endif

