/*
Duo Zhao
*/

#include "particleSystem.h"
#include "showSystem.h"
#include "input.h"
#include "physics.h"

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 30;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int shear=0, bend=0, structural=1, pause=1, viewingMode=1, saveScreenToFile=0;

struct world jello;
struct particleSystem cps;


int windowWidth, windowHeight;

void initialCPS(){
	
  cps.integrator = 0;
  cps.dt = 0.001; // timestep, e.g.. 0.001
  cps.n = 5; // display only every nth timepoint
  cps.mass = 10; // mass of each of particle
  cps.resolution = 0; // resolution for the 3d grid specifying the external force field; value of 0 means that there is no force field
  
	cps.p[0].x=0;cps.p[0].y=12;cps.p[0].z=0;
	cps.p[1].x=0;cps.p[1].y=8;cps.p[1].z=0;
	cps.p[2].x=0;cps.p[2].y=4;cps.p[2].z=0;
	cps.p[3].x=0;cps.p[3].y=0;cps.p[3].z=0;
	cps.p[4].x=4;cps.p[4].y=0;cps.p[4].z=0;
	cps.p[5].x=8;cps.p[5].y=0;cps.p[5].z=0;
	cps.p[6].x=12;cps.p[6].y=0;cps.p[6].z=0;
  for(int i=0; i<N_particle; i++){
	cps.v[i].x=3;cps.v[i].y=3;cps.v[i].z=0;   
  }

}

void myinit()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(90.0,1.0,0.01,1000.0);

  // set background color to grey
  glClearColor(0.5, 0.5, 0.5, 0.0);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  return; 
}

void reshape(int w, int h) 
{
  // Prevent a divide by zero, when h is zero.
  // You can't make a window of zero height.
  if(h == 0)
    h = 1;

  glViewport(0, 0, w, h);

  // Reset the coordinate system before modifying
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Set the perspective
  double aspectRatio = 1.0 * w / h;
  gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); 

  windowWidth = w;
  windowHeight = h;

  glutPostRedisplay();
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // camera parameters are Phi, Theta, R
  
  gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
	        0.0,0.0,0.0, 0.0,0.0,1.0);

  /* Lighting */
  /* You are encouraged to change lighting parameters or make improvements/modifications
     to the lighting model . 
     This way, you will personalize your assignment and your assignment will stick out. 
  */

  // global ambient light
  GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };
  
  // light 's ambient, diffuse, specular
  GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd0[] = { 0.0, 1.0, 0.0, 1.0 };
  GLfloat lKs0[] = { 0.0, 1.0, 0.0, 1.0 };

  GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd1[] = { 1.0, 0.0, 0.0, 1.0 };
  GLfloat lKs1[] = { 1.0, 0.0, 0.0, 1.0 };

  GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd2[] = { 0.0, 0.0, 1.0, 1.0 };
  GLfloat lKs2[] = { 0.0, 0.0, 1.0, 1.0 };

  GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd3[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKs3[] = { 0.0, 1.0, 1.0, 1.0 };

  GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd4[] = { 1.0, 0.0, 1.0, 1.0 };
  GLfloat lKs4[] = { 1.0, 0.0, 1.0, 1.0 };

  GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd5[] = { 1.0, 1.0, 0.0, 1.0 };
  GLfloat lKs5[] = { 1.0, 1.0, 0.0, 1.0 };

  GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd6[] = { 0.3, 0.7, 0.3, 1.0 };
  GLfloat lKs6[] = { 0.3, 0.7, 0.3, 1.0 };

  GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd7[] = { 0.7, 0.3, 0.7, 1.0 };
  GLfloat lKs7[] = { 0.7, 0.3, 0.7, 1.0 };
  

  // light positions and directions
  GLfloat lP0[] = { 30.0, 30.0, 30.0, 1.0 };
  GLfloat lP1[] = { 30.0, 30.0, -30.0, 1.0 };
  GLfloat lP2[] = { -30.0, 30.0, 30.0, 1.0 };
  GLfloat lP3[] = { -30.0, 30.0, -30.0, 1.0 };
  GLfloat lP4[] = { 30.0, -30.0, 30.0, 1.0 };
  GLfloat lP5[] = { 30.0, -30.0, -30.0, 1.0 };
  GLfloat lP6[] = { -30.0, -30.0, 30.0, 1.0 };
  GLfloat lP7[] = { -30.0, -30.0, -30.0, 1.0 };
  
  // jelly material color

  GLfloat mKa[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat mKd[] = { 0.5, 0.5, 0.5, 1.0 };
  GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 };

  /* set up lighting */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  // set up cube color
  glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
  glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
  glMaterialf(GL_FRONT, GL_SHININESS, 120);
    
  // macro to set up light i
  #define LIGHTSETUP(i)\
  glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
  glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
  glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
  glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
  glEnable(GL_LIGHT##i)
  
  LIGHTSETUP (0);
  LIGHTSETUP (1);
  LIGHTSETUP (2);
  LIGHTSETUP (3);
  LIGHTSETUP (4);
  LIGHTSETUP (5);
  LIGHTSETUP (6);
  LIGHTSETUP (7);

  // enable lighting
  glEnable(GL_LIGHTING);    
  glEnable(GL_DEPTH_TEST);

   //show the constrained Particle System
  showCPS(&cps);

  glDisable(GL_LIGHTING);
  glutSwapBuffers();
}

void doIdle()
{
  char s[20]="picxxxx.ppm";  
  // save screen to file
  s[3] = 48 + (sprite / 1000);
  s[4] = 48 + (sprite % 1000) / 100;
  s[5] = 48 + (sprite % 100 ) / 10;
  s[6] = 48 + sprite % 10;

  if (saveScreenToFile==1&&pause == 0)
  {
    saveScreenshot(windowWidth, windowHeight, s);
    saveScreenToFile=1; // save only once, change this if you want continuos image generation (i.e. animation)
    sprite++;
  }

  if (sprite >= 300) // allow only 300 snapshots
  {
    exit(0);	
  }

  if (pause == 0)
  {
    // insert code which appropriately performs one step of the cube simulation:
	  
	  Euler(&cps);
	 
  }//end of pause


  glutPostRedisplay();
}

int main (int argc, char ** argv)
{
  initialCPS();
  Euler(&cps);
  glutInit(&argc,argv);
  
  /* double buffered window, use depth testing, 640x480 */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  
  windowWidth = 640;
  windowHeight = 480;
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitWindowPosition (0,0);
  glutCreateWindow ("Constrained Particle System");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mouseMotionDrag);

  /* callback for window size changes */
  glutReshapeFunc(reshape);

  /* callback for mouse movement */
  glutPassiveMotionFunc(mouseMotion);

  /* callback for mouse button changes */
  glutMouseFunc(mouseButton);

  /* register for keyboard events */
  glutKeyboardFunc(keyboardFunc);

  /* do initialization */
  myinit();

  /* forever sink in the black hole */
  glutMainLoop();

  return(0);
}

