/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "particleSystem.h"
#include "physics.h"

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */
  //go through each mass point
  int i,j,k;
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
		{
		  //compute F_total, initialize it
		  struct point F_total;
		  F_total.x=0;F_total.y=0;F_total.z=0;
		  struct point L; //Vector pointing from B to A
		  struct point F_temp; //Each Force
		  struct point V; //Va-Vb
		  double L_length;
		  double scale;
		  int flag=0;
		  
		  //for each mass point, first check if it's inside box
		  //printf("f   %f\n", abs(6.5));
		  //printf("d   %d\n", abs(6));

		  double distanceToBoundingBox=0;
		  if(fabs(jello->p[i][j][k].x)>2.0){
			  flag=1;
			  distanceToBoundingBox=fabs(jello->p[i][j][k].x)-2.0;
			  if(jello->p[i][j][k].x>0){
				  L.x=1; L.y=0; L.z=0;
			  }else{
				  L.x=-1; L.y=0; L.z=0;
			  }
		  }else
		  if(fabs(jello->p[i][j][k].y)>2.0){
			  flag=1;
			  distanceToBoundingBox=fabs(jello->p[i][j][k].y)-2.0;
			  if(jello->p[i][j][k].y>0){
				  L.x=0; L.y=1; L.z=0;
			  }else{
				  L.x=0; L.y=-1; L.z=0;
			  }
		  }else
		  if(fabs(jello->p[i][j][k].z)>2.0){
			  flag=1;
			  distanceToBoundingBox=fabs(jello->p[i][j][k].z)-2.0;
			  if(jello->p[i][j][k].z>0){
				  L.x=0; L.y=0; L.z=1;
			  }else{
				  L.x=0; L.y=0; L.z=-1;
			  }
		  }
		  if(flag==1){
		  //compute L
		  pMULTIPLY(L,distanceToBoundingBox*DISTANCE_FACTOR,L);
		  //compute collision spring force
		  scale = -jello->kCollision;
		  pMULTIPLY(L, scale, F_temp);
		  pSUM(F_total, F_temp, F_total);
		  //compute collision damping force
		  V.x=jello->v[i][j][k].x; 
		  V.y=jello->v[i][j][k].y; 
		  V.z=jello->v[i][j][k].z;
		  scale = -jello->dCollision*(V.x*L.x+V.y*L.y+V.z*L.z)/distanceToBoundingBox/distanceToBoundingBox;
		  pMULTIPLY(L, scale, F_temp);
		  pSUM(F_total, F_temp, F_total);
		  }
		  //inclined plane
		  if(jello->incPlanePresent==1){
			  double lineEquation = jello->a*jello->p[i][j][k].x+jello->b*jello->p[i][j][k].y+jello->c*jello->p[i][j][k].z+jello->d;
			  double distanceToPlane=0;
			  if(lineEquation<0){
				  distanceToPlane=lineEquation/sqrt(pow(jello->a,2)+pow(jello->b,2)+pow(jello->c,2));
				  L.x=jello->a;  L.y=jello->b;  L.z=jello->c;
				  double length;
				  pNORMALIZE(L);
				  pMULTIPLY(L,distanceToPlane*DISTANCE_FACTOR,L);
				  pMULTIPLY(L, -jello->kCollision, F_temp);
				  pSUM(F_total, F_temp, F_total);
			  }
		  }
		  

		  //loop through 27 nearby mass points
		  for(int ii=-1; ii<=1; ii++)
			  for(int jj=-1; jj<=1; jj++)
				  for(int kk=-1; kk<=1; kk++)
					{
					  //check if the current point is valid
					  if(i+ii<0||i+ii>7||j+jj<0||j+jj>7||k+kk<0||k+kk>7){
						  continue;
					  }
					  //find the distance, determine the spring length
					  if((abs(ii)+abs(jj)+abs(kk))==0){
						  //current point
						  continue;
					  }
					  if((abs(ii)+abs(jj)+abs(kk))==1){
						  //distance=R_structural
						  pDIFFERENCE(jello->p[i][j][k],jello->p[i+ii][j+jj][k+kk],L);
						  L_length = sqrt((L).x * (L).x + (L).y * (L).y + (L).z * (L).z);
						  scale = -jello->kElastic*(L_length-R_structural)/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
						  //Damping force
						  pDIFFERENCE(jello->v[i][j][k],jello->v[i+ii][j+jj][k+kk],V);
						  scale = -jello->dElastic*(V.x*L.x+V.y*L.y+V.z*L.z)/L_length/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
						  //when sum equals one, there are bend springs
						  //bend spring F_hook
						  int itemp,jtemp,ktemp;
						  itemp=jtemp=ktemp=0;
						  if(ii==1){
							 itemp=2;
						  }else if(ii==-1){
							 itemp=-2;
						  }else if(jj==1){
							 jtemp=2;
						  }else if(jj==-1){
							 jtemp=-2;
						  }else if(kk==1){
							 ktemp=2;
						  }else if(kk==-1){
							 ktemp=-2;
						  }
						  if(!(i+itemp<0||i+itemp>7||j+jtemp<0||j+jtemp>7||k+ktemp<0||k+ktemp>7)){
						  //Hook force
						  pDIFFERENCE(jello->p[i][j][k],jello->p[i+itemp][j+jtemp][k+ktemp],L);
						  L_length = sqrt((L).x * (L).x + (L).y * (L).y + (L).z * (L).z);
						  scale = -jello->kElastic*(L_length-R_bend)/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
						  
						  //Damping force
						  pDIFFERENCE(jello->v[i][j][k],jello->v[i+itemp][j+jtemp][k+ktemp],V);
						  scale = -jello->dElastic*(V.x*L.x+V.y*L.y+V.z*L.z)/L_length/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
					    }
						  
					  }
					  if((abs(ii)+abs(jj)+abs(kk))==2){
						  //distance=R_shear_short
						  pDIFFERENCE(jello->p[i][j][k],jello->p[i+ii][j+jj][k+kk],L);
						  L_length = sqrt((L).x * (L).x + (L).y * (L).y + (L).z * (L).z);
						  scale = -jello->kElastic*(L_length-R_shear_short)/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
						  //Damping force
						  pDIFFERENCE(jello->v[i][j][k],jello->v[i+ii][j+jj][k+kk],V);
						  scale = -jello->dElastic*(V.x*L.x+V.y*L.y+V.z*L.z)/L_length/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
					  }
					  if((abs(ii)+abs(jj)+abs(kk))==3){
						  //distance=R_shear_long
						  pDIFFERENCE(jello->p[i][j][k],jello->p[i+ii][j+jj][k+kk],L);
						  L_length = sqrt((L).x * (L).x + (L).y * (L).y + (L).z * (L).z);
						  scale = -jello->kElastic*(L_length-R_shear_long)/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
						  //Damping force
						  pDIFFERENCE(jello->v[i][j][k],jello->v[i+ii][j+jj][k+kk],V);
						  scale = -jello->dElastic*(V.x*L.x+V.y*L.y+V.z*L.z)/L_length/L_length;
						  pMULTIPLY(L, scale, F_temp);
						  pSUM(F_total, F_temp, F_total);
					  }  
				    }
					
			//apply forcefield
			//find nearby 8 forcefield sample and do interpolation
				  int resolution=jello->resolution;
				  if(flag==0){ //current point is inside boundbox, x y z from -2 to 2
					  double forceFieldUnitLength = 4.0/(resolution-1);
					  double x,y,z;
					  x=jello->p[i][j][k].x+2.0;
					  y=jello->p[i][j][k].y+2.0;
					  z=jello->p[i][j][k].z+2.0; //now x y z from 0 to 4
					  int indexX, indexY, indexZ;
					  indexX=(x-fmod(x,forceFieldUnitLength))/forceFieldUnitLength;
					  indexY=(y-fmod(y,forceFieldUnitLength))/forceFieldUnitLength;
					  indexZ=(z-fmod(z,forceFieldUnitLength))/forceFieldUnitLength;	

					  if(indexX>=0&&indexX<=resolution-2&&indexY>=0&&indexY<=resolution-2&&indexZ>=0&&indexZ<=resolution-2){
						  double average1, average2, average3, average4, average12, average34, average1234x, average1234y, average1234z;
						  
						  //interpolate x
						  average1 = jello->forceField[(indexX+0)*resolution*resolution+(indexY+0)*resolution + indexZ+0].x*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+0)*resolution*resolution+(indexY+0)*resolution + indexZ+1].x*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average2 = jello->forceField[(indexX+1)*resolution*resolution+(indexY+0)*resolution + indexZ+0].x*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+1)*resolution*resolution+(indexY+0)*resolution + indexZ+1].x*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average3 = jello->forceField[(indexX+0)*resolution*resolution+(indexY+1)*resolution + indexZ+0].x*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+0)*resolution*resolution+(indexY+1)*resolution + indexZ+1].x*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average4 = jello->forceField[(indexX+1)*resolution*resolution+(indexY+1)*resolution + indexZ+0].x*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+1)*resolution*resolution+(indexY+1)*resolution + indexZ+1].x*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average12 = average1*((forceFieldUnitLength-fmod(x,forceFieldUnitLength))/forceFieldUnitLength)
											+average2*(fmod(x,forceFieldUnitLength)/forceFieldUnitLength);
						  average34 = average3*((forceFieldUnitLength-fmod(x,forceFieldUnitLength))/forceFieldUnitLength)
											+average4*(fmod(x,forceFieldUnitLength)/forceFieldUnitLength);
						  
						   F_temp.x = average12*((forceFieldUnitLength-fmod(y,forceFieldUnitLength))/forceFieldUnitLength)
									    +average34*(fmod(y,forceFieldUnitLength)/forceFieldUnitLength);

						  //interpolate y
						  average1 = jello->forceField[(indexX+0)*resolution*resolution+(indexY+0)*resolution + indexZ+0].y*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+0)*resolution*resolution+(indexY+0)*resolution + indexZ+1].y*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average2 = jello->forceField[(indexX+1)*resolution*resolution+(indexY+0)*resolution + indexZ+0].y*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+1)*resolution*resolution+(indexY+0)*resolution + indexZ+1].y*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average3 = jello->forceField[(indexX+0)*resolution*resolution+(indexY+1)*resolution + indexZ+0].y*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+0)*resolution*resolution+(indexY+1)*resolution + indexZ+1].y*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average4 = jello->forceField[(indexX+1)*resolution*resolution+(indexY+1)*resolution + indexZ+0].y*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+1)*resolution*resolution+(indexY+1)*resolution + indexZ+1].y*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average12 = average1*((forceFieldUnitLength-fmod(x,forceFieldUnitLength))/forceFieldUnitLength)
											+average2*(fmod(x,forceFieldUnitLength)/forceFieldUnitLength);
						  average34 = average3*((forceFieldUnitLength-fmod(x,forceFieldUnitLength))/forceFieldUnitLength)
											+average4*(fmod(x,forceFieldUnitLength)/forceFieldUnitLength);
						  
						   F_temp.y = average12*((forceFieldUnitLength-fmod(y,forceFieldUnitLength))/forceFieldUnitLength)
									    +average34*(fmod(y,forceFieldUnitLength)/forceFieldUnitLength);

						  //interpolate z
						  average1 = jello->forceField[(indexX+0)*resolution*resolution+(indexY+0)*resolution + indexZ+0].z*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+0)*resolution*resolution+(indexY+0)*resolution + indexZ+1].z*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average2 = jello->forceField[(indexX+1)*resolution*resolution+(indexY+0)*resolution + indexZ+0].z*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+1)*resolution*resolution+(indexY+0)*resolution + indexZ+1].z*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average3 = jello->forceField[(indexX+0)*resolution*resolution+(indexY+1)*resolution + indexZ+0].z*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+0)*resolution*resolution+(indexY+1)*resolution + indexZ+1].z*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average4 = jello->forceField[(indexX+1)*resolution*resolution+(indexY+1)*resolution + indexZ+0].z*((forceFieldUnitLength-fmod(z,forceFieldUnitLength))/forceFieldUnitLength)
							        +jello->forceField[(indexX+1)*resolution*resolution+(indexY+1)*resolution + indexZ+1].z*(fmod(z,forceFieldUnitLength)/forceFieldUnitLength);
						  average12 = average1*((forceFieldUnitLength-fmod(x,forceFieldUnitLength))/forceFieldUnitLength)
											+average2*(fmod(x,forceFieldUnitLength)/forceFieldUnitLength);
						  average34 = average3*((forceFieldUnitLength-fmod(x,forceFieldUnitLength))/forceFieldUnitLength)
											+average4*(fmod(x,forceFieldUnitLength)/forceFieldUnitLength);
						  
						   F_temp.z = average12*((forceFieldUnitLength-fmod(y,forceFieldUnitLength))/forceFieldUnitLength)
									    +average34*(fmod(y,forceFieldUnitLength)/forceFieldUnitLength);

						   //add force field to total force
						   pSUM(F_total, F_temp, F_total);
					  }
				  }
				  
			//update a
		  a[i][j][k].x=F_total.x/jello->mass;
		  a[i][j][k].y=F_total.y/jello->mass;
		  a[i][j][k].z=F_total.z/jello->mass;
		}
	  
}

void computeAcceleration(struct particleSystem * cps, struct point a[N_particle]){
	for(int i=0; i<N_particle; i++){
		a[i].x=0; a[i].y=0;
	}

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

