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
#include <iostream>
#include <fstream> // for file handling

#include <QDebug> // USED FOR DEBUGGING

// const defining number of points in one circle of a sphere
#define CIRCLE_POINTS 100 

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
			p.x = cosPhi * cosTheta;
			p.y = cosPhi * sinTheta;
			p.z = sinPhi;
			spherePoints.push_back(p); 
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
	
	// load modelview for translation and scaling
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix(); // push matrix to apply object specific transforms
	glTranslatef(centre.x, centre.y, centre.z);
	glScalef(radius, radius, radius);	
			
	// draw the strips of our sphere 
	glBegin(GL_QUAD_STRIP);
	// create CIRCLE_POINTS-1 quad primitive latitude strips
	for(int i = 0; i < (CIRCLE_POINTS - 1); i++)  {
		/* creates quads between one sphere and sphere above has to include
		the first point in sphere at start and at end hence why the for loop
		does CIRCLE_POINTS+1 iterations */
		for(int j = 0; j < (CIRCLE_POINTS + 1); j++)  {
			// get points
			Point lower = spherePoints.at((i * CIRCLE_POINTS) + (j % CIRCLE_POINTS));
			Point upper = spherePoints.at(((i+1) * CIRCLE_POINTS) + (j % CIRCLE_POINTS));
			// render vertices
			glVertex3d(lower.x, lower.y, lower.z);
			glVertex3d(upper.x, upper.y, upper.z);
		}
	}			
	glEnd();
	
	// fills the hole at the bottom of the sphere
	glBegin(GL_POLYGON);
	for(int i = 0; i < CIRCLE_POINTS; i++) {
		Point p = spherePoints.at(i); // get point
		glVertex3d(p.x, p.y, p.z); // render vertex
	}
	glEnd();
	
	// fills the hole at the top of the sphere
	glBegin(GL_POLYGON);
	for(int i = 0; i < CIRCLE_POINTS; i++){
		Point p = spherePoints.at(CIRCLE_POINTS * (CIRCLE_POINTS - 1) + i); // get point
		glVertex3d(p.x, p.y, p.z); // render vertex
	}
	glEnd();
	
	glPopMatrix(); // return to normal modelview matrix
}

/**
   writes centre and radius information from file to vectors
	   
   @param objData file containing object raw data 
   @return returns 1 if there is a major error, else 0
*/
int vis::setData(char* objData) {
	
	double x, y, z, rad; // will temporarily hold file inputs
	int i = 0; // counter for error output
	
	ifstream dataFile(objData);
	if (dataFile.is_open()) { // opens file if location is valid
	
		/* data is read in this way so if a value is missed from the file it can be detected.
		   the x is read first as an eof check, then the error counter is incremented so we can
		   compare this later with the size of our vectors. If a value is missing, there will 
		   be a discrepency between the error counter and vector size */
		while(dataFile >> x) {
			i++; // adds to the increment counter, used for error testing
			if(dataFile >> y >> z >> rad) { // reads lines of file
				if(rad > 0.0) { // only radius over 0 valid
					// write to vector
					objCentre.push_back((Point){x, y, z}); 
					objRadius.push_back(rad);
					cout << x << " " << y << " " << z << " " << rad << "\n";
				
				} else { // radius error
					cout << "ERROR: line " << i << " of file is invalid and will be ignored!\n";
					i -= 1; // takes one off counter to ignore line
				}
			} 
		}
		
		dataFile.close();		
		
		// check if file hasincorrect number of doubles specified
		if((int)objRadius.size() < i) {
			cout << "ERROR: centre and radius not specified for all objects!\n";
			return 1; // exits qt program
		}
		
		// check file wasn't empty
		if(i == 0) {
			cout << "ERROR: file was empty!\n";
			return 1; // exits qt program
		}
		
	} else { // invalid file specified
		cout << "ERROR: file specified is invalid!\n";
		return 1; // exits qt program
	}
	
	return 0; // succesful return
}	

/**
   initialises environment for OpenGL rendering when instance called
*/
void vis::initializeGL() {

	// initialise sphere points
	spherePoints = makeSpherePoints(CIRCLE_POINTS);
	glClearColor(0.0,0.0,0.0,0.0); // black background
	glMatrixMode(GL_PROJECTION); // projection mode to set clipping plane
	glLoadIdentity();
	glFrustum(-2.0,2.0,-1.0,1.0,2.0,1000.0); // sets the clipping plane
	glTranslatef(0.0, 0.0, -10.0); // moves camera back to view scene
	glMatrixMode(GL_MODELVIEW); // initialise modelview matrix
	glLoadIdentity();
	scaleFactor = 1.0; // initialise the zoom factor variable
	
	glEnable(GL_DEPTH_TEST); // allows for depth comparison when renderin
	
	// setting up lighting
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0); // creating one light, light0
	GLfloat fPosition[4] = {0.0, 1.0, 0.0, 1.0}; // light position
	glLightfv(GL_LIGHT0, GL_POSITION, fPosition); // setting position
	// now to specify ambient, diffuse, and specular intensities, all white
	GLfloat fiAmbient[4] = {0.2, 0.2, 0.2, 1.0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, fiAmbient); 
	GLfloat fiDiffuse[4] = {0.8, 0.8, 0.8, 1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, fiDiffuse);
	GLfloat fiSpecular[4] = {0.0, 0.0, 0.0, 0.0};
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
	
	// draw objects from vectors
	for(int i = 0; i < (int)objCentre.size(); i++)
		drawSphere(objCentre.at(i), objRadius.at(i));
	
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
	
	glMatrixMode(GL_PROJECTION); // change to modelview matrix
	if(deg > 0.0) {// wheel scrolled forward
		glScalef(pow(1.1, deg), pow(1.1, deg), pow(1.1, deg)); 
		scaleFactor = scaleFactor*pow(1.1, deg); // update scale factor
	}
	else {// wheel scrolled backward
		glScalef(pow(0.9, -deg), pow(0.9, -deg), pow(0.9, -deg));
		scaleFactor = scaleFactor*pow(0.9, -deg); // update scale factor
	}
	
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
	
	glMatrixMode(GL_PROJECTION); // change to modelview matrix
	if(key == Qt::Key_Equal) { // "+" button pressed
		glScalef(1.1f, 1.1f, 1.1f); // zoom in 
		scaleFactor = scaleFactor*1.1; // update scale factor
	}
	else if(key == Qt::Key_Minus) {// "-" button pressed
		glScalef(0.9f, 0.9f, 0.9f); // zoom out
		scaleFactor = scaleFactor*0.9; // update scale factor
	}
		
	event->accept(); // accepts the event
	updateGL(); // redraw to screen
}

/**
   interaction handling for when a button on the mouse is pressed,  
   used for camera translation and rotation, and picking. 
   
   @param event information about the mouse button press
*/
void vis::mousePressEvent(QMouseEvent *event) { 
	startPos = event->pos(); // records position that mouse was clicked
}

void vis::mouseMoveEvent(QMouseEvent *event) {
	
	// calculate change in x and y from start point to current mouse pos
	float xPos = (float)(event->x() - startPos.x());
	float yPos = (float)(event->y() - startPos.y());
	
	// translate camera case, translates based on current zoom scale
	if(event->buttons() == Qt::RightButton) {
	
		glMatrixMode(GL_MODELVIEW); // load modelview matrix to be translated
		glTranslatef((1.0/scaleFactor)*(xPos/68.0), (1.0/scaleFactor)*(-yPos/68.0), 0.0); // translate view based on zoom scale
		
	} else if (event->buttons() == Qt::LeftButton) {
	
		glMatrixMode(GL_MODELVIEW); // load projection matrix to be translated
		// rotate x direction movement about y axis
		glRotatef(xPos, 0.0, 1.0, 0.0); 
		// NEEED TO CHANGE THIS!
		// rotate y direction whatever vector is perpendicular to screen
		glRotatef(yPos, 1.0, 0.0, 0.0);
		
	}
	
	event->accept(); // accepts event
	updateGL(); // redraw screen
	startPos = event->pos(); // update the start position
	
}
