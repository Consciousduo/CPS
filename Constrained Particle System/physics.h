

#ifndef _PHYSICS_H_
#define _PHYSICS_H_
#define DISTANCE_FACTOR 5.0


void computeAcceleration(struct particleSystem * cps, struct point a[N_particle]);

void printMatrixC(int row, int column, double matrix[][14]);
void printMatrixT(int row, int column, double matrix[][9]);
void printMatrixA(int row, int column, double matrix[][23]);

void Euler(struct particleSystem * cps);
#endif

