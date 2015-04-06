//show constrined particle system

#include "particleSystem.h"
#include "showSystem.h"

void showCPS(struct particleSystem * cps){
	//show particles
	for(int i=0; i<N_particle; i++){ 
		glTranslatef(cps->p[i].x,cps->p[i].y,cps->p[i].z);
		glutSolidSphere( 1.0, 50.0, 50.0);
		glTranslatef(-cps->p[i].x,-cps->p[i].y,-cps->p[i].z);
	}
	//glDisable(GL_LIGHTING);
	glLineWidth(2.5); 
	glColor3f(0.0, 1.0, 0.0);
	//show connections
	for(int i=0; i<N_particle-1; i++){
		glBegin(GL_LINES);
		glVertex3f(cps->p[i].x, cps->p[i].y, cps->p[i].z);
		glVertex3f(cps->p[i+1].x, cps->p[i+1].y, cps->p[i+1].z);
		glEnd();
	}
	//show circle
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINE_LOOP);
	for(int i =0; i <= 300; i++){
	double angle = 2 * pi * i / 300;
	double x = cos(angle)*12;
	double y = sin(angle)*12;
	glVertex2d(x,y);
	}
	glEnd();

}

