/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "particleSystem.h"
#include "physics.h"

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */

void computeAcceleration(struct particleSystem * cps, struct point a[N_particle]){
	for(int i=0; i<N_particle; i++){
		a[i].x=0; a[i].y=0;
	}

/*
  double a_data[] = { 0.18, 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, 0.58,
                      0.14, 0.30, 0.97, 0.66,
                      0.51, 0.13, 0.19, 0.85 };

  double b_data[] = { 1.0, 2.0, 3.0, 4.0 };

  gsl_matrix_view m 
    = gsl_matrix_view_array (a_data, 4, 4);

  gsl_vector_view b
    = gsl_vector_view_array (b_data, 4);

  gsl_vector *x = gsl_vector_alloc (4);
  
  int s;

  gsl_permutation * p = gsl_permutation_alloc (4);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

  printf ("x = \n");
  gsl_vector_fprintf (stdout, x, "%g");

  gsl_permutation_free (p);
  gsl_vector_free (x);
*/
	//construct mass matrix mass_data[row][column]
	double mass_data[N_particle*2][N_particle*2];
	for(int i=0; i<N_particle*2; i++){
		for(int j=0; j<N_particle*2; j++){
			if(i!=j)
			mass_data[i][j]=0;
			else
			mass_data[i][j]=cps->mass;
		}
	}
	//construct ¦¤C matrix mass_data[row][column]
	double ¦¤C_data[N_particle+1][N_particle*2];
	//construct ¦¤CT matrix mass_data[row][column]
	double ¦¤CT_data[N_particle*2][N_particle+1];

	for(int i=0; i<N_particle*2; i++){
		¦¤C_data[0][i] = 0;
		¦¤C_data[1][i] = 0;
	}
	¦¤C_data[0][0] = 1; ¦¤C_data[1][1] = 1; 

	for(int i=2; i<N_particle+1; i++){ //2 ~ row
		for(int j=0; j<N_particle*2; j++){ //all column
			int k = i-2;
			if(j == (i-2)*2){
				¦¤C_data[i][j] = 2*(cps->p[k].x-cps->p[k+1].x);
			}
			else if(j == (i-2)*2+1){
				¦¤C_data[i][j] = 2*(cps->p[k].y-cps->p[k+1].y);
			}
			else if(j == (i-2)*2+2){
				¦¤C_data[i][j] = 2*(cps->p[k+1].x-cps->p[k].x);
			}
			else if(j == (i-2)*2+3){
				¦¤C_data[i][j] = 2*(cps->p[k+1].y-cps->p[k].y);
			}
			else{
				¦¤C_data[i][j]=0;
			}
		}
	}

	for(int i=0; i<N_particle*2; i++){
		for(int j=0; j<N_particle+1; j++){
			¦¤CT_data[i][j] = ¦¤C_data[j][i];
		}
	}
	//printMatrixC(N_particle+1, N_particle*2, ¦¤C_data);
	//printMatrixT(N_particle*2, N_particle+1, ¦¤CT_data);
	//construct the final matrix
	double A_data[3*N_particle+1][3*N_particle+1];
	for(int i=0; i<3*N_particle+1; i++){
		for(int j=0; j<3*N_particle+1; j++){
			if(i<2*N_particle && j<2*N_particle){
				A_data[i][j] = mass_data[i][j];
			}
			else if(i<2*N_particle && j>=2*N_particle){
				A_data[i][j] = ¦¤CT_data[i][j-2*N_particle];
			}
			else if(i>=2*N_particle && j<2*N_particle){
				A_data[i][j] = ¦¤C_data[i-2*N_particle][j];
			}else{
				A_data[i][j]=0;
			}

		}
	}
	//A_data finished
	//printMatrixA(3*N_particle+1, 3*N_particle+1, A_data);



}


/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct particleSystem * cps)
{
  int i,j,k;
  point a[N_particle];
  computeAcceleration(cps, a);
  for(int i=0; i<N_particle; i++){
		cps->v[i].x += cps->dt * a[i].x;
        cps->v[i].y += cps->dt * a[i].y;
		cps->p[i].x += cps->dt * cps->v[i].x;
        cps->p[i].y += cps->dt * cps->v[i].y;
  }
  
 
}


void printMatrixC(int row, int column, double matrix[][14])
{
	for(int i=0; i<row; i++){
		for(int j=0; j<column; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
}

void printMatrixT(int row, int column, double matrix[][8])
{
	for(int i=0; i<row; i++){
		for(int j=0; j<column; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
}

void printMatrixA(int row, int column, double matrix[][22])
{
	for(int i=0; i<row; i++){
		for(int j=0; j<column; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
}
