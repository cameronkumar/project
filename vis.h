/**
   Class declaration for the visualisation part of the non-Euclidean 3D
   visualisation software. An instance of this will appear in the GUI 
   
   @author Cameron Kumar
   @version 1.0 5/3/16
*/

#include <QGLWidget> // allows us to define a Qt OpenGL widget
#include <vector>

// simplified calling functions from the std class
using namespace std;

// header guard, stops repeated inclusion of this header by the compiler
#ifndef __VIS_H__
#define __VIS_H__

// defines a struct to hold triples of doubles for use as coordinate points
struct Point {
	
	double x, y, z;	
	
};

// defines a structure to hold quads of floats for use as RGBA values
struct RGBA {

	float R, G, B, A;

};

// defines a structure that contains an id and distance from plane value
// used for sorting translucent objects 
struct idDist {

	int id;
	double dist;

};

// defines a structure that contains an id, overlap distance, direction vector and length
// used for sorting storing intersection details
struct idOverVecLen {
	
	int id;
	double over;
	Point vec;
	double len;	
	
};

// structure to store data about a particular intersection, including centre of intersection details
// and the ids of the two spheres involved. Used for drawing intersections
struct intDraw {

	Point cen;
	double rad;
	int id1;
	int id2;

};

// defining our visualisation class as subclass of QGLWidget
class vis: public QGLWidget {

	public:
	
	/**
	   class constructor
	   
	   @param parent parent class, inherit from there 
	*/
	vis(QWidget *parent);
	
	/**
	   creates the points of a sphere object, centre at origin, radius 1
	   
	   @param nPoints number of points per circle in the sphere
	   @return vector of double triplets defining the point coords in 3d 
	*/
	vector<Point> makeSpherePoints(int nPoints);
	
	/**
	   draws a sphere with specified centre and radius
	   
	   @param centre coordinate point of centre location
	   @param radius radius of sphere 
	*/
	void drawSphere(Point centre, double radius);
	
	/**
	   writes centre and radius information from file to vectors
	   
	   @param objData file containing object raw data 
	   @return returns 1 if there is a major error, else 0
	*/
	int setData(char* objData);
	
	/**
	   populate the colours vector with a selection of RGB values
	*/
	void initColours();
	
	/**
	   Recursive function to be used for merge sorting of translucent spheres
	   
	   @param vec vector to sort (structure of ints and doubles)
	   @return sorted vector (structure of ints and doubles)
	*/
	vector<idDist> mergeSort(vector<idDist>);
	
	
	/**
	   orders the spheres based on opacity and which is closest to the camera
	   
	   @return returns a vector of integers representing the order of furthest to nearest spheres
	*/
	vector<int> sphereOrder();
	
	/**
	   Changes the RGB colour of an object defined by id to the RGB provided
	   
	   @param id identifier of object whose colour we want to change
	   @param rgb colour we want to change the object to
	*/
	void changeColour(int id, Point rgb);
	
	/**
	   Changes the transparency value of object to alpha value provided
	   
	   @param id identifier of object whose transparency we want to change
	   @param alpha alpha value we want to change object to
	*/
	void changeTransparency(int id, double alpha);
	
	/**
	   Returns human readable string detailing all intersections for an object
	   
	   @param id index of object to calculate intersections for
	   @param inter list of intersecting objects
	   @return human readable string containing intersection details
	*/
	string getIntersectionString(int id, vector<idOverVecLen> inter); 
	
	/**
	   Returns a string detailing all intersections and tangents for specified object
	   
	   @param id identifier of object we will calculate intersections for
	   @return string in human readable form of intersection details
	*/
	string intersectsWith(int id);
	
	/**
	   Draws a circle given a specified centre and radius
   
   	   @param cen centre of circle to draw
   	   @param rad radius of circle to draw
	*/
	void drawCircle(intDraw circ);
	
	/** 
	   Draws the intersections saved in the global variable coi
	*/
	void drawIntersections();
	
	private:
	
	// holds the coordinate points for a standard sphere
	vector<Point> spherePoints; 
	
	// holds the data of all geometric objects centres
	vector<Point> objCentre;
	
	// holds the radius data of geometric objects
	vector<double> objRadius;
	
	// holds the RGB colour values of each object (default red), and transparency value A (default 1.0)
	vector<RGBA> objColour;
	
	// holds the RGB colour options for objects
	vector<Point> colour;
	
	// holds the starting position when the mouse is clicked
	QPoint startPos;
	
	// holds the current zoom scale of the camera, used to scale translation
	double scaleFactor;
	
	// holds current rotation about the yaw axis, used for calculating transforms
	double yRot;
	
	// holds the current rotation about the pitch axis, used to prevent over rotation
	double pRot;
	
	// holds information about the centre of circle of intersections for intersections to be drawn
	vector<intDraw> coi;

	protected:
	
	/**
	   initialises environment for OpenGL rendering when instance called
	*/
	void initializeGL();
	
	/**
	   changes size of viewport when widget resized
	  
	   @param w the new width of the widget
	   @param h the new height of the widget
	*/
	void resizeGL(int w, int h);
	
	/**
	   draw a new frame
	*/ 
	void paintGL();
	
	/**
	   interaction handling for the mouses's scroll wheel, used for 
	   camera zoom
	   
	   @param event information about the mouse button click
	*/
	void wheelEvent(QWheelEvent *event);
	
	/**
	   interaction handling for when a button on the keyboard is pressed,  
	   used for camera zooming
	   
	   @param event information about the key pressed
	*/
	void keyPressEvent(QKeyEvent *event);
	
	/**
	   interaction handling for when a button on the mouse is pressed,  
	   used for camera translation and rotation, and picking
	   
	   @param event information about the mouse button press
	*/
	void mousePressEvent(QMouseEvent *event);
	
	/**
	   interaction handling for when mouse movements take place, used for
	   camera translations and rotations
	   
	   @param event information about the mouse button press
	*/
	void mouseMoveEvent(QMouseEvent *event);
	
}; // end of class

#endif // ends our header guard
