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
void vis::drawSphere(Point centre, double radius) {
	
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
			// render vertices with normals, as the sphere points are for a sphere
			// with radius 1 and centre (0,0,0), normal vector is just each sphere point
			// interpretted as a vector
			glNormal3d(lower.x, lower.y, lower.z); 
			glVertex3d(lower.x, lower.y, lower.z);
			glNormal3d(upper.x, upper.y, upper.z);
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
					// set default colour and transparency values
					objColour.push_back((RGBA){1.0, 0.0, 0.0, 1.0});
				
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
   populate the colours vector with a selection of RGB values
*/
void vis::initColours() {
	
	colour.push_back((Point){1.0, 0.0, 0.0}); // red (default)
	colour.push_back((Point){0.0, 1.0, 0.0}); // green 
	colour.push_back((Point){0.0, 0.0, 1.0}); // blue 
	colour.push_back((Point){0.0, 1.0, 1.0}); // cyan	
	colour.push_back((Point){1.0, 0.0, 1.0}); // pink
	colour.push_back((Point){1.0, 1.0, 0.0}); // yellow
	colour.push_back((Point){1.0, 1.0, 1.0}); // white
	
}

/**
   Recursive function to be used for merge sorting of translucent spheres
   
   @param vec vector to sort (structure of ints and doubles)
   @return sorted vector (structure of ints and doubles)
*/
vector<idDist> vis::mergeSort(vector<idDist> vec) {
	
	// terminates when only 1 element in list
	if(vec.size() == 1)
		return vec;
	else {
	
		// determine middle position and left and right vectors
		vector<idDist>::iterator mid = vec.begin() + (vec.size()/2);
		vector<idDist> l(vec.begin(), mid);
		vector<idDist> r(mid, vec.end());
		
		// merge sort each component (recursive step)
		l = mergeSort(l);
		r = mergeSort(r);
		
		// sort the left and right vectors into an order
		vector<idDist> sorted; // will hold sorted values
		int iL = 0, iR = 0; // increment variables for left and right loops
		
		// adding sorted objects to the sorted vector by comparison until end of one vector reached
		while(iL < (int)l.size() && iR < (int)r.size()) {
			if(l[iL].dist < r[iR].dist) {
				sorted.push_back(l[iL]);
				iL++;
			} else {
				sorted.push_back(r[iR]);
				iR++;
			}
		}
		
		// now adding remaining contents of vectors to sorted list
		while(iL < (int)l.size()) {
			sorted.push_back(l[iL]);
			iL++;
		}
		while(iR < (int)r.size()) {
			sorted.push_back(r[iR]);
			iR++;
		}
		
		return sorted; // returning the completed sorted vector
	}
	
}

/**
   orders the spheres based on opacity and which is closest to the camera
	   
   @return returns a vector of integers representing the order of furthest to nearest spheres
*/
vector<int> vis::sphereOrder() {
	
	vector<idDist> translucentid; // will hold id and distance of translucent objects
	vector<int> opaqueid; // will hold id of opaque objects
	
	// need to seperate opaque from translucent objects first before we sort translucents
	for(int i = 0; i < (int)objCentre.size(); i++) {
		if(objColour.at(i).A == 1.0)
			opaqueid.push_back(i);
		else
			translucentid.push_back((idDist){i, 0.0});
	}
	
	// initialise vector that will store final order list, opaque objects rendered first
	vector<int> order = opaqueid; 
	
	// if there are translucent objects to be ordered then that is done now
	if((int)translucentid.size() > 0) {
		// to sort translucent objects, first calculate the distance from camera
		// we will work out order based on point-plane distance from the cameras
		// parallel plane at origin, need to work out plane normal first
		Point norm = {(sin((-yRot*M_PI)/180.0))*cos((pRot*M_PI)/180.0), // x
			      sin((pRot*M_PI)/180.0), // y
			      cos((-yRot*M_PI)/180.0)*cos((pRot*M_PI)/180.0)}; // z
		// we also need to calculate the normal vectors magnitude for the equation
		double normMag = sqrt(pow(norm.x, 2) + pow(norm.y, 2) + pow(norm.z, 2)); 
	
		for(int i = 0; i < (int)translucentid.size(); i++) {
	
			Point cen = objCentre.at(translucentid.at(i).id); // get centre
			double rad = objRadius.at(translucentid.at(i).id); // get radius
		
			// calculate signed distance from plane
			double d = (norm.x*cen.x + norm.y*cen.y + norm.z*cen.z)/normMag;
			translucentid.at(i).dist = d + rad; // update vector
	
		}
	
		// now need to sort the vector, for this we use mergesort recursion
		translucentid = mergeSort(translucentid);
	
		// finally recollate ordered list of ids to draw
		for(int i = 0; i < (int)translucentid.size(); i++) 
			order.push_back(translucentid.at(i).id);
	} 

	return order; // returning ordered list

}

/**
   Returns human readable string detailing all intersections for an object
	   
   @param id index of object to calculate intersections for
   @param inter list of intersecting objects
   @return human readable string containing intersection details
*/
string vis::getIntersectionString(int id, vector<idDist> inter) {

}

/**
   Returns a string detailing all intersections and tangents for specified object
	   
   @param id identifier of object we will calculate intersections for
   @return string in human readable form of intersection details
*/
string vis::intersectsWith(int id) {
	
	// get the object's centre and radius
	Point c = objCentre.at(id); 
	double r = objRadius.at(id);
	// create a vector to store details of which objects intersect
	vector<idDist> inter;
	
	// this for loop calculates which intersections occur
	for(int i = 0; i < (int)objCentre.size(); i++) {
		if(i != id) { // don't want to compare with self
		
			// get centre and radius for current i
			Point iC = objCentre.at(i);
			double iR = objRadius.at(i);
			
			// calculate length between spheres and compare to sum of radii
			double length = sqrt(pow(c.x-iC.x,2) + pow(c.y-iC.y,2) + pow(c.z-iC.z,2));
			double overlap = (r+iR) - length; 
			
			// determining if an intersection occurs, if so added to vector
			if(overlap >= 0.0)
				inter.push_back((idDist){i, overlap});
					
		}	
	}
	
	// now to create string to be returned to user
	return(getIntersectionString(int id, vector<idDist> inter)); 
}

/**
   Changes the RGB colour of an object defined by id to the RGB provided
   
   @param id identifier of object whose colour we want to change
   @param rgb colour we want to change the object to
*/
void vis::changeColour(int id, Point rgb) {
	objColour.at(id).R = rgb.x;
	objColour.at(id).G = rgb.y;
	objColour.at(id).B = rgb.z;
}

/**
   Changes the transparency value of object to alpha value provided
	   
   @param id identifier of object whose transparency we want to change
   @param alpha alpha value we want to change object to
*/
void vis::changeTransparency(int id, double alpha) {
	objColour.at(id).A = alpha;
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
	pRot = yRot = 0.0; // initialise rotation variables
	initColours(); // initialise colour vector
	changeColour(2, colour.at(2));
	changeColour(3, colour.at(1));
	changeTransparency(2, 0.4);
	changeTransparency(3, 0.7);

	
	glEnable(GL_DEPTH_TEST); // allows for depth comparison when renderin
	
	// setting up lighting
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0); // creating one light, light0
	// now to specify ambient and diffuse intensities, all white
	GLfloat fiAmbient[4] = {0.2, 0.2, 0.2, 1.0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, fiAmbient); 
	GLfloat fiDiffuse[4] = {0.8, 0.8, 0.8, 1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, fiDiffuse);
	// set global ambient lighting to illuminate scene more
	GLfloat global_ambient[4] = {0.5, 0.5, 0.5, 1.0};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	
	// enable blending, this will be used for translucent objects
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
}

/**
   changes size of viewport when widget resized
 
   @param w the new width of the widget
   @param h the new height of the widget
*/
void vis::resizeGL(int w, int h) {
	
	glViewport(0, 0, w, h);
	
	float ratio = (float)w/(float)h;
	
	// now change the frustum to reflect the width and height;
	glMatrixMode(GL_PROJECTION); // projection mode to set clipping plane
	glLoadIdentity();
	glFrustum(-ratio,ratio,-1.0,1.0,2.0,1000.0); // sets the clipping plane
	glTranslatef(0.0, 0.0, -10.0); // moves camera back to view scene
	glScalef(scaleFactor, scaleFactor, scaleFactor);

}

/**
   draws a new frame
*/ 
void vis::paintGL() {

	// create the render order
	vector<int> renderOrder = sphereOrder();

	// clearing screen first
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// initialise standard surface properties and light position for all objects
	GLfloat fPosition[4] = {1.0, 1.0, 1.0, 0.0}; // light position
	glLightfv(GL_LIGHT0, GL_POSITION, fPosition); // setting light position
	
	// draw objects from vectors, specify surface properties from vector
	for(int i = 0; i < (int)renderOrder.size(); i++) {
		// get colour from vector and set properties
		RGBA col = objColour.at(renderOrder.at(i));
		glMaterialfv(GL_FRONT, GL_AMBIENT, (GLfloat*)&col);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, (GLfloat*)&col);
		// draw the sphere
		drawSphere(objCentre.at(renderOrder.at(i)), objRadius.at(renderOrder.at(i)));
	}
	
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
	
	// translate camera case, translates based on current zoom scale and rotation about y
	// we translate in 2 dimensions, up and down in y direction, and left and right in axis
	// parallel to screen. This axis direction needs to be calculated from y rotation
	if(event->buttons() == Qt::RightButton) {
	
		glMatrixMode(GL_MODELVIEW); // load modelview matrix to be translated
		glTranslatef(cos((yRot*M_PI)/180.0)*(1.0/scaleFactor)*(xPos/68.0), // x direction translation
		             (1.0/scaleFactor)*(-yPos/68.0), // y direction translation
		             sin((yRot*M_PI)/180.0)*(1.0/scaleFactor)*(xPos/68.0)); // z direction translation
		
	} else if (event->buttons() == Qt::LeftButton) {
	
		glMatrixMode(GL_MODELVIEW); // load modelview matrix to be translated
		// rotate x direction movement about yaw axis
		glRotatef(xPos, 0.0, 1.0, 0.0); 
		// we want to keep track of the rotation about the yaw axis, so update var
		yRot += xPos;
		
		// rotate about pitch axis if rotation doesnt exceed limits
		if(-90.0<pRot+yPos && pRot+yPos<90.0) {
			glRotatef(yPos, cos((yRot*M_PI)/180.0), 0.0, sin((yRot*M_PI)/180.0));
			// we want to keep track of pitch rotation so it doesnt go over 180
			pRot += yPos;
		}
		
	}
	
	event->accept(); // accepts event
	updateGL(); // redraw screen
	startPos = event->pos(); // update the start position
	
}
