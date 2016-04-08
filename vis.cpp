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
#include <string>
#include <sstream>
#include <iomanip>

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
					objColour.push_back((RGBA){1.0, 0.0, 0.0, 0.4});
				
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
   Returns human readable string detailing all intersections for an object
	   
   @param id index of object to calculate intersections for
   @param inter list of intersecting objects
   @return human readable string containing intersection details
*/
string vis::getIntersectionString(int id, vector<idOverVecLen> inter) {
	
	stringstream sInter; // will hold human readable intersection details
	sInter << setprecision(4); // setting precision for doubles
	// get the centre coordinates and radius for calculations
	Point cen = objCentre.at(id); 
	double rad = objRadius.at(id);
	
	if((int)inter.size() == 0) // case of no intersections
		sInter << "Sphere " << id << " does not intersect any other objects!\n";
	else { 
		// loop through all intersections and append string with information
		for(int i = 0; i < (int)inter.size(); i++) {
			
			// read intersection details from vector for calculations
			idOverVecLen currentInt = inter.at(i);
		
			if(inter.at(i).over == 0.0) { // tangent case
			
				// calculating tangent point
				Point pTangent = (Point){cen.x+rad*currentInt.vec.x, // x
						 	 cen.y+rad*currentInt.vec.y, // y
						 	 cen.z+rad*currentInt.vec.z};// z
				// write to intersection string
				sInter << "Sphere " << id << " is tangent to sphere " << currentInt.id << " at point ("
				       << pTangent.x << ", " << pTangent.y << ", " << pTangent.z << ")\n";
				// add tangent to intersection vector for drawing, relatively small radius
				coi.push_back((intDraw){pTangent, currentInt.vec, objRadius.at(i)*0.01, id, i});
				       
				
			} else { 
				
				// get radius of intersecting sphere
				double iRad = objRadius.at(currentInt.id);
				// now calculate the distance from the centre of sphere i to coi 
				double coiDist = pow(currentInt.len, 2) - (pow(iRad,2) - pow(rad,2));
				coiDist = coiDist/(2*currentInt.len);
				// also calculate the radius using the distance calculation as this is a part of the equation
				double coiRad = sqrt(pow(rad,2) - pow(coiDist,2));
				// use this to calculate the centre of the circle
				Point coiCen = (Point){cen.x + coiDist*currentInt.vec.x,
					       	       cen.y + coiDist*currentInt.vec.y,
					       	       cen.z + coiDist*currentInt.vec.z};
				
				// case of tangency between one sphere inside another
				if(rad == iRad + currentInt.len || rad == iRad - currentInt.len) { 
					
					// write tangency string
					sInter << "Sphere " << id << " is tangent to sphere " << currentInt.id << " at point ("
					       << coiCen.x << ", " << coiCen.y << ", " << coiCen.z << ")\n";
					// add to vector
					coi.push_back((intDraw){coiCen, currentInt.vec,  objRadius.at(i)*0.01, id, i}); 		
					
					// write to intersection string
					
				} else if(rad < iRad + currentInt.len) { // intersection case
				
					// write to intersection string
					sInter << "Sphere " << id << " intersects sphere " << currentInt.id 
				       	       << " with circle of intersection located about (" << coiCen.x << ", " 
				               << coiCen.y << ", " << coiCen.z << ") with radius " << coiRad << "\n";	
					// add intersection to intersection vector for drawing
					coi.push_back((intDraw){coiCen, currentInt.vec, coiRad, id, i});
				
				}
				
				// nothing done if sphere completely contained in another sphere
				 			
			}
		}	
	}
	return(sInter.str()); // return string stream string
}

/**
   Returns a string detailing all intersections and tangents for specified object
	   
   @param id identifier of object we will calculate intersections for
   @return string in human readable form of intersection details
*/
string vis::intersectsWith(int id) {
	
	// error check to ensure valid id specified
	if(id < 0 || id >= (int)objCentre.size())
		return("Error! ID specified for intersection check out of range!\n");
	else {
		// get the object's centre and radius
		Point c = objCentre.at(id); 
		double r = objRadius.at(id);
		// create a vector to store details of which objects intersect
		vector<idOverVecLen> inter;
	
		// this for loop calculates which intersections occur
		for(int i = 0; i < (int)objCentre.size(); i++) {
			if(i != id) { // don't want to compare with self
		
				// get centre and radius for current i
				Point iC = objCentre.at(i);
				double iR = objRadius.at(i);
			
				// calculate length between spheres, normalize direction vector and compare to sum of radii
				Point cenVec = (Point){iC.x - c.x, iC.y - c.y, iC.z - c.z};
				double cenVecLength = sqrt(pow(cenVec.x,2) + pow(cenVec.y,2) + pow(cenVec.z,2));
				cenVec = (Point){cenVec.x/cenVecLength, cenVec.y/cenVecLength, cenVec.z/cenVecLength};
				double overlap = (r+iR) - cenVecLength; 
			
				// determining if an intersection occurs, if so added to vector
				if(overlap >= 0.0)
					inter.push_back((idOverVecLen){i, overlap, cenVec, cenVecLength});
					
			}	
		}
	
		// now to create string to be returned to user
		return(getIntersectionString(id, inter)); 
	}
}

/**
   Draws a circle given a specified circle of intersection
   
   @param circ information about circle to draw
*/
void vis::drawCircle(intDraw circ) {
	
	// get the direction vector as this will be orthogonal vector of coi
	Point dirVec = circ.vec;	
	
	// set up the translation matrix to move the circle to correct location
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix(); 
	glTranslatef(circ.cen.x, circ.cen.y, circ.cen.z); // translate to centre point
	glScalef(circ.rad, circ.rad, circ.rad);	// scale to radius size
	/* current orthogonal vector (0, 0, 1), we need to pitch & yaw rotate to dirVec,
	   this can be done using vector equations*/
	glRotatef((180*acos(dirVec.z)/M_PI), -dirVec.y, dirVec.x, 0.0); // cross product 
	
	// standard circle of with radius 1 about (0, 0, 0) on y x-y plane
	// as translation matrix has been set up above this will be drawn at position & orientation
	glBegin(GL_LINE_LOOP);
	for(int i = 0; i < CIRCLE_POINTS; i++) {		
		
		double angle = ((double)i/(double)CIRCLE_POINTS)*2*M_PI; // calc current angle
			
		// calc current point on polygon
		Point p = (Point){sin(angle), cos(angle), 0.0};
		
		// add point to GL_POLYGON
		glNormal3d(p.x, p.y, p.z);
		glVertex3d(p.x, p.y, p.z);
		
	}
	glEnd();
	
	glPopMatrix(); // return to normal modelview matrix
	
}

/** 
   Draws the intersections saved in the global variable coi
*/
void vis::drawIntersections() {

	// for each intersection
	for(int i = 0; i < (int)coi.size(); i++) {
		// get details of circle of intersection from vector
		intDraw iCoi = coi.at(i);
		// draw circle
		drawCircle(iCoi);
	}
	
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
	trans = 0; // initialise translation flag
	select = 0; // initialise selection flag
	
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
	
	for(int i = 0; i < (int)objCentre.size(); i++) {
		objColour.at(i).R = colour.at(i%6).x;
		objColour.at(i).G = colour.at(i%6).y;
		objColour.at(i).B = colour.at(i%6).z;
	}
	
	cout << intersectsWith(0);
	cout << intersectsWith(1);
	cout << intersectsWith(2);
	cout << intersectsWith(3);
	cout << intersectsWith(4);
	cout << intersectsWith(5);
	cout << intersectsWith(6);
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
	
	// drawing intersections first, draw them solid white
	RGBA white = (RGBA){1.0, 1.0, 1.0, 1.0};
	glMaterialfv(GL_FRONT, GL_AMBIENT, (GLfloat*)&white);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, (GLfloat*)&white);
	drawIntersections();
	
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
	
		trans = 1; // set flag to indicate a camera translation is taking place
	
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

/**
   interaction handling for when mouse button released, used for picking
	   
   @param event information about the mouse button released
*/
void vis::mouseReleaseEvent(QMouseEvent *event) {

	if(event->buttons() == Qt::RightButton) { 
		
		if(trans==1) // case of translation currently taking place
			trans = 0; // reset translation flag to 0 as it has ended
		else { // case of picking
			select = 1; // set selection mode flag
			updateGL(); // render in selection mode
		}	
		
	}
	
}
