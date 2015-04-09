/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "particleSystem.h"
#include "physics.h"
#include <gsl/gsl_linalg.h>


/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */

void computeAcceleration(struct particleSystem * cps, struct point a[N_particle]){
	

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
	double ¦¤C_data[N_particle+2][N_particle*2];
	//construct ¦¤CT matrix mass_data[row][column]
	double ¦¤CT_data[N_particle*2][N_particle+2];

	for(int i=0; i<N_particle*2; i++){
		¦¤C_data[0][i] = 0;
		¦¤C_data[1][i] = 0;
	}
	¦¤C_data[0][0] = 1; ¦¤C_data[1][1] = 1; 

	for(int i=2; i<N_particle+1; i++){ //row
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
	//first n+1 row finished
	for(int i=0; i<N_particle*2; i++){
		¦¤C_data[N_particle+1][i]=0;
		if(i==(N_particle*2-2))
			¦¤C_data[N_particle+1][i]=2*cps->p[N_particle-1].x;
		if(i==(N_particle*2-1))
			¦¤C_data[N_particle+1][i]=2*cps->p[N_particle-1].y;
	}
	//all row finished

	for(int i=0; i<N_particle*2; i++){
		for(int j=0; j<N_particle+2; j++){
			¦¤CT_data[i][j] = ¦¤C_data[j][i];
		}
	}
	//printMatrixC(N_particle+2, N_particle*2, ¦¤C_data);
	//printMatrixT(N_particle*2, N_particle+2, ¦¤CT_data);
	//construct the final matrix
	double A_data[3*N_particle+2][3*N_particle+2];
	for(int i=0; i<3*N_particle+2; i++){
		for(int j=0; j<3*N_particle+2; j++){
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
	//printMatrixA(3*N_particle+2, 3*N_particle+2, A_data);

	//construct dc'/dq matrix, we will reuse ¦¤C_data since they are the same size
	for(int i=0; i<N_particle+2; i++){
		for(int j=0; j<N_particle*2; j++){
			¦¤C_data[i][j] = 0;
		}
	}
	for(int i=2; i<N_particle+1; i++){
		for(int j=0; j<N_particle*2; j++){
			int k = i-2;
			if(j == (i-2)*2){
				¦¤C_data[i][j] =-2*(cps->v[k+1].x-cps->v[k].x);
			}
			else if(j == (i-2)*2+1){
				¦¤C_data[i][j] =-2*(cps->v[k+1].y-cps->v[k].y);
			}
			else if(j == (i-2)*2+2){
				¦¤C_data[i][j] = 2*(cps->v[k+1].x-cps->v[k].x);
			}
			else if(j == (i-2)*2+3){
				¦¤C_data[i][j] = 2*(cps->v[k+1].y-cps->v[k].y);
			}
			else{
				¦¤C_data[i][j]=0;
			}
		}
	}	
	¦¤C_data[N_particle+1][N_particle*2-2]=2*cps->v[N_particle-1].x;
	¦¤C_data[N_particle+1][N_particle*2-1]=2*cps->v[N_particle-1].y;
	//end of construction of dc'/dq matrix
	//printMatrixC(N_particle+2, N_particle*2, ¦¤C_data);

	//construct q'
	double qp[N_particle*2];
	for(int i=0; i<N_particle; i++){
		qp[2*i]  = cps->v[i].x;
		qp[2*i+1]= cps->v[i].x;
	}

	//construct dc'/dq * q'
	double result[N_particle+2];
	for(int i=0; i<N_particle+2; i++){
		double temp_sum=0;
		for(int j=0; j<N_particle*2; j++){
			temp_sum+=¦¤C_data[i][j]*qp[j];
		}
		result[i]=-temp_sum;//add the negative sign here
	}
	/*
	for(int i=0; i<N_particle+2; i++){
		printf("%f\n", result[i]);
	}
	*/
	//Substract Baumgarte Stabilization from result
	//Construct C
	double C[N_particle+2];
	C[0] = cps->p[0].x;
	C[1] = cps->p[0].y-12;
	C[N_particle+1] = cps->p[N_particle-1].x*cps->p[N_particle-1].x + cps->p[N_particle-1].y*cps->p[N_particle-1].y - 144;
	for(int i=2; i<N_particle+1; i++){
		int k = i-2;
		double t1 = cps->p[k+1].x-cps->p[k].x;
		double t2 = cps->p[k+1].y-cps->p[k].y;
		C[i] = t1*t1+t2*t2-16.1;//here we can adjust the distance constrain of each sphere 4.2*4.2
	} // now we get C
	/*
	for(int i=0 ;i<N_particle+2; i++){
		printf("%f\n", C[i]);
	}
	*/
	//Construc C'
	double Cp[N_particle+2];
	Cp[0] = cps->v[0].x;
	Cp[1] = cps->v[0].y;
	Cp[N_particle+1] = 2*cps->p[N_particle-1].x*cps->v[N_particle-1].x + 2*cps->p[N_particle-1].y*cps->v[N_particle-1].y;
	for(int i=2; i<N_particle+1; i++){
		int k = i-2;
		Cp[i] = 2*(cps->p[k+1].x-cps->p[k].x)*(cps->v[k+1].x-cps->v[k].x) + 2*(cps->p[k+1].y-cps->p[k].y)*(cps->v[k+1].y-cps->v[k].y);
	} 
	/*
	for(int i=0 ;i<N_particle+2; i++){
		printf("%f\n", Cp[i]);
	}*/
	//Substract Baumgarte Stabilization from result
	for(int i=0; i<N_particle+2; i++){
		double b = cps->b;
		result[i] = result[i] - b*b*C[i] - 2*b*Cp[i];
	}

	//finally Construct b vector; //here we add gravity and any other forces
	double B_data[3*N_particle+2];
	for(int i=0; i<N_particle; i++){
		B_data[i*2] = 0-cps->vDamping*cps->v[i].x;           //x
		B_data[i*2+1] = -9.8 - cps->vDamping*cps->v[i].y;;   //y
	}
	for(int i=2*N_particle; i<N_particle*3+2; i++){
		B_data[i] = result[i-2*N_particle];
	}
	//B_data is ready
	/*
	for(int i=0; i<3*N_particle+2; i++){
		printf("%f\n", B_data[i]);
	}*/
	double a_data[(3*N_particle+2)*(3*N_particle+2)];
	for(int i=0; i<(3*N_particle+2); i++){
		for(int j=0; j<(3*N_particle+2); j++){
			a_data[i*(3*N_particle+2)+j] = A_data[i][j];
		}
	}

	for(int i=0; i<N_particle; i++){
		a[i].x=0; 
		a[i].y=0;
	}
    
	  
  gsl_matrix_view leftSide = gsl_matrix_view_array (a_data, 3*N_particle+2, 3*N_particle+2);
  gsl_vector_view rightSide = gsl_vector_view_array (B_data, 3*N_particle+2);
  gsl_vector *x = gsl_vector_alloc (3*N_particle+2); 
  int s;
  gsl_permutation * p = gsl_permutation_alloc (3*N_particle+2);
  gsl_linalg_LU_decomp (&leftSide.matrix, p, &s);
  gsl_linalg_LU_solve (&leftSide.matrix, p, &rightSide.vector, x);
  
  for(int i=0; i<N_particle; i++){
	    a[i].x=x->data[i*2]; 
		a[i].y=x->data[i*2+1];
	}

  gsl_permutation_free (p);
  gsl_vector_free (x);
  
}


/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct particleSystem * cps)
{
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

	printf("\n");printf("\n");
}

void printMatrixT(int row, int column, double matrix[][9])
{
	for(int i=0; i<row; i++){
		for(int j=0; j<column; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}

	printf("\n");printf("\n");
}

void printMatrixA(int row, int column, double matrix[][23])
{
	for(int i=0; i<row; i++){
		for(int j=0; j<column; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}

	printf("\n");printf("\n");
}

