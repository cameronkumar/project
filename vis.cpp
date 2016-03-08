/**
   Function definitions for the visualisation class of the non-Euclidean 3D
   visualisation software. All member function code is defined here 
   
   @author Cameron Kumar
   @version 1.0 5/3/16
*/

#include "vis.h"
#include <QMouseEvent>
#include <QWheelEvent>
#include <QKeyEvent>
#include <cmath> // for power & other math funcs

#include <QDebug> // USED FOR DEBUGGING

// const defining number of points in one circle of a sphere
#define CIRCLE_POINTS 100; 

/**
   class constructor
   
   @param parent parent class, inherit from there 
*/
vis::vis(QWidget *parent): QGLWidget(parent) {} // simple constuctor

/**
   creates the points of a sphere object, centre at origin, radius 1
   
   @param nPoints number of points per circle in the sphere
   @return vector of Point structure defining the point coords in 3d 
*/
vector<Point> vis::makeSpherePoints(int nPoints) {
	
	vector<Point> spherePoints; // temp vector for points
	
	// loops to create the sphere points as circles
	for(int i = 0; i < nPoints; i++) {
		
		// calculations for latitude location of circle of points
		float phi = (-(M_PI) / 2.0) + (((float)i / (float)nPoints) * M_PI);
		float cosPhi = cos(phi); // used for calculation
		float sinPhi = sin(phi); // used for calculation
		
		for(int j = 0; j < nPoints; j++) {
			
			Point p; // temp point
			
			// calculate longtitude value of each point in circle
			float theta = -(M_PI) + (((float)j / (float)nPoints) * 2.0 * M_PI);
			float cosTheta = cos(theta); // used for calculation
			float sinTheta = sin(theta); // used for calculation
			
			// now to add coord to vector
			p.x =  radius * cosPhi * cosTheta;
			p.y = radius * cosPhi * sinTheta;
			p.z = radius * sinPhi;
			spherePoints.pushback(p); 
		}
	}
	
	return spherePoints;
	
}

/**
   draws a sphere with specified centre and radius
   
   @param centre coordinate point of centre location
   @param radius radius of sphere 
*/
void vis::drawSphere(Point centre, double radius){
	
	// draw the strips of our sphere 
	glBegin(GL_QUAD_STRIP);

	// create CIRCLE_POINTS-1 quad primitive latitude strips
	for(int i = 0; i < CIRCLE_POINTS - 1; i++)  {
		
		/* creates quads between one sphere and sphere above has to include
		the first point in sphere at start and at end hence why the for loop
		does CIRCLE_POINTS+1 iterations */
		for(int j = 0; j < CIRCLE_POINTS + 1; j++)  {
		
			glVertex3d(spherePoints.at((i * CIRCLE_POINTS) + (j % CIRCLE_POINTS)));
			glVertex3d(spherePoints.at(((i+1) * CIRCLE_POINTS) + (j % CIRCLE_POINTS)));
			
		}
		
	}
				
	glEnd();
	
	// fills the hole at the bottom of the sphere
	glBegin(GL_POLYGON);
	
	for(int i = 0; i < CIRCLE_POINTS; i++)
		glVertex3d(spherePoints.at(i));
		
	glEnd();
	
	// fills the hole at the top of the sphere
	glBegin(GL_POLYGON);
	
	for(int i = 0; i < CIRCLE_POINTS; i++)
		glVertex3d(spherePoints.at((CIRCLE_POINTS * (CIRCLE_POINTS - 1) + i));
		
	glEnd();
	
}


/**
   initialises environment for OpenGL rendering when instance called
*/
void vis::initializeGL() {

	// initialise sphere points
	spherePoints = makeSpherePoints(CIRCLE_POINTS);
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
   draws a new frame
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
	
	Point p = {0.0, 0.0, 0.0};
	
	drawSphere(p, 1.0);
	
	// draw frame and render to screen
	glFinish();

}

/**
   interaction handling for the mouses's scroll wheel, used for 
   camera zoom
   
   @param event information about the mouse button click
*/
void vis::wheelEvent(QWheelEvent *event) {

	// deg +120 for one roll forward, -120 for one roll back, read this here
	float deg = (float)event->delta();
	deg = deg/120.0; // deg now number of scrolls in + or - dir'n
	
	glMatrixMode(GL_MODELVIEW); // change to modelview matrix

	if(deg > 0.0) // wheel scrolled forward
		glScalef(pow(1.1, deg), pow(1.1, deg), pow(1.1, deg)); 
	else // wheel scrolled backward
		glScalef(pow(0.9, -deg), pow(0.9, -deg), pow(0.9, -deg));
	
	event->accept(); // accepts the event
	updateGL(); // redraw to screen
}

/**
   interaction handling for when a button on the keyboard is pressed,  
   used for camera zooming
   
   @param event information about the key pressed
*/
void vis::keyPressEvent(QKeyEvent *event) {
	
	int key = event->key(); // get the integer value of key pressed
	
	glMatrixMode(GL_MODELVIEW); // change to modelview matrix
	
	if(key == Qt::Key_Equal) // "+" button pressed
		glScalef(1.1f, 1.1f, 1.1f); // zoom in 
	else if(key == Qt::Key_Minus) // "-" button pressed
		glScalef(0.9f, 0.9f, 0.9f); // zoom out
		
	event->accept(); // accepts the event
	updateGL(); // redraw to screen
}
