/**
   Function definitions for the visualisation class of the non-Euclidean 3D
   visualisation software. All member function code is defined here 
   
   @author Cameron Kumar
   @version 1.0 5/3/16
*/

#include "vis.h"
#include <QMouseEvent>
#include <QKeyEvent>

/**
   class constructor
   
   @param parent parent class, inherit from there 
*/
vis::vis(QWidget *parent): QGLWidget(parent) {} // simple constuctor

/**
   initialises environment for OpenGL rendering when instance called
*/
void vis::initializeGL() {

	glClearColor(0.0,0.0,0.0,0.0); // black background
	glOrtho(-1.0,1.0,-1.0,1.0,-1.0,1.0); // sets the clipping plane
	
	glEnable(GL_DEPTH_TEST); // allows for depth comparison when renderin
	
	// setting up lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0); // creating one light, light0
	GLfloat fPosition[4] = {0.0, 0.0, 0.0, 1.0}; // light position
	glLightfv(GL_LIGHT0, GL_POSITION, fPosition); // setting position
	// now to specify ambient, diffuse, and specular intensities, all white
	GLfloat fiAmbient[4] = {0.4, 0.4, 0.4, 1.0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, fiAmbient); 
	GLfloat fiDiffuse[4] = {0.6, 0.6, 0.6, 1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, fiDiffuse);
	GLfloat fiSpecular[4] = {1.0, 1.0, 1.0, 1.0};
	glLightfv(GL_LIGHT0, GL_SPECULAR, fiSpecular);
	
}

/**
   changes size of viewport when widget resized
 
   @param w the new width of the widget
   @param h the new height of the widget
*/
void vis::resizeGL(int w, int h) {

	// resize the viewport to that of the new widget dimensions
	glViewport(0,0,w,h);
	
}

/**
   draw a new frame
*/ 
void vis::paintGL() {

	// clearing screen first
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// initialise default surface properties for all objects, red colour
	GLfloat fRed[4] = {1.0, 0.0, 0.0, 1.0};
	GLfloat fWhite[4] = {1.0, 1.0, 1.0, 1.0};
	glMaterialfv(GL_FRONT, GL_AMBIENT, fRed); 
	glMaterialfv(GL_FRONT, GL_DIFFUSE, fRed); 
	glMaterialfv(GL_FRONT, GL_SPECULAR, fWhite); 
	glMaterialf(GL_FRONT, GL_SHININESS, 128.0);
	// temporarily drawing a primitive before spheres implemented 
	glBegin(GL_QUADS);
		glVertex3f(0.5f, 0.5f, 0.0f);
		glVertex3f(-0.5f, 0.5f, 0.0f);
		glVertex3f(-0.5f, -0.5f, 0.0f);
		glVertex3f(0.5f, -0.5f, 0.0f);
	glEnd();
	// draw frame and render to screen
	glFinish();

}

/**
   interaction handling for when a mouse button is pressed, used for 
   camera movement
   
   @param event information about the mouse button click
*/
void vis::mousePressEvent(QMouseEvent *event) {}

/**
   interaction handling for when a button on the keyboard is pressed,  
   used for camera zooming
   
   @param event information about the key pressed
*/
void vis::keyPressEvent(QKeyEvent *event) {}
