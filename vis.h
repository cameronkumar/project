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

// defines a structure to store generation and key data of objects, will be used 
// with a vector to store information for all objects
struct genKey {

	int gen;
	char key;
	
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

// structure to store data about a particular intersection, including centre of intersection details,
// vector between the two spheres, and the ids of the two spheres involved. Used for drawing intersections
struct intDraw {

	Point cen, vec;
	double rad;
	int id1, id2;

};

// stores 2 integers, used when pair cannot be for non-compiletime constants
struct intPair {
	
	int x, y;
	
};

// stores data about the parameters of slideshow mode specified by the user
struct slideParam {

	char key;
	int trans, group;

};

// defining our visualisation class as subclass of QGLWidget
class vis: public QGLWidget {

	// Q_OBJECT macro
	Q_OBJECT

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
	   create points of a cube about the origin with side length 1
	   cube points added front to back, top to bottom, left to right
	   
	   @return vector containing points of cube
	*/
	vector<Point> makeCubePoints();
	
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
	   populate the colours vector with a selection of RGB values and creates colour icons
	*/
	void initColours();
	
	/**
	   Determines whether a value is already present within a list or not
	   
	   @param val value to be searched for
	   @param list list to be searched
	   @return 1 if true, else 0
	*/
	int valPresent(char val, vector<char> list);
	
	/**
	   Finds all unique key values of objects
	   
	   @return vector containing all unique key values
	*/
	vector<char> getKeyList();
	
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
	void changeColour(int id, RGBA rgb);
	
	/**
	   Changes the transparency value of object to alpha value provided
	   
	   @param id identifier of object whose transparency we want to change
	   @param alpha alpha value we want to change object to
	*/
	void changeTransparency(int id, double alpha);
	
	/**
	   Calculates the centre, radius, and orthogonal vector for an intersection
	   Stores detail in global vector
	   
	   @param id index of object to calculate intersections for
	   @param inter list of intersecting objects
	*/
	void calculateIntersection(int id, vector<idOverVecLen> inter); 
	
	/**
	   Identifies which objects intersect specified object
	   
	   @param id identifier of object we will calculate intersections for
	   @return number of objects this object intersects with
	*/
	int intersectsWith(int id);
	
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
	
	/**
	   Renders for picking
	   
	   @return integer id of picked object
	*/
	int getPicked();
	
	/**
	   Draws cube side with face length 1 about origin
	   
	   @param p array of int containing each point on cube face
	*/
	void drawCubeFace(int p[4]);
	
	/**
	   Renders a wireframe cube around picked object
	   
	   @param id of sphere that cube should be drawn around
	*/
	void selectionCube(int id);
	
	/**
	   Creates the context menu for a right click pick
	*/
	void createContextMenu();
	
	/**
	   Creates the dialog pop-up for scroll mode
	*/
	void createScrollDialog();
	
	/**
	   Scroll mode command to move selection to next object
	   
	   @param delta the rotation amount of the scroll wheel
	*/
	void scroll(float delta);
	
	/**
	   Function to print specific intersection between two objects and add interstion
	   to draw list
	   
	   @param a index of first sphere
	   @param b index of second sphere
	*/
	void handleIntersection(int a, int b);
	
	/** 
	   resets all variables back to normal when slideshow mode ends
	*/
	void stopSlideshow();
	
	/**
	   Function to set up drawing of next slide to screen in slideshow mode
	*/
	void createSlide();

	
	private:
	
	// define number of points in one circle of a sphere
	int circPoints;
	
	// holds the coordinate points for a standard sphere
	vector<Point> spherePoints; 
	
	// holds the coordinate points for standard cube
	vector<Point> cubePoints;
	
	// holds the data of all geometric objects centres
	vector<Point> objCentre;
	
	// holds the radius data of geometric objects
	vector<double> objRadius;
	
	// holds the RGB colour values of each object (default red), and transparency value A (default 1.0)
	vector<RGBA> objColour;
	
	// holds the generation and key value for each object
	vector<genKey> objGenKey;
	
	// holds the RGB colour options for objects
	vector<RGBA> colour;
	
	// vector to store all unqiue key valeus
	vector<char> key;
	
	// holds the QIcon types for each colour to be used in combobox
	vector<QIcon*> colourIcon;
	
	// holds the starting position when the mouse is clicked
	QPoint startPos;
	
	// holds the current zoom scale of the camera, used to scale translation
	double scaleFactor;
	
	// holds current rotation about the yaw axis, used for calculating transforms
	double yRot;
	
	// holds the current rotation about the pitch axis, used to prevent over rotation
	double pRot;
	
	// holds the position that a particular objects cois begin at in the coi vector
	vector<int> coiBegin;
	
	// holds information about circles of intersection
	vector<intDraw> coi;
	
	// holds information about circles of intersection to be drawn
	vector<intDraw> coiDraw; 
	
	// flag to indicate whether a camera translation is currently taking place
	int trans;
	
	// flag to indicate whether we wish to render for colour picking
	int colourPicking;
	
	// stores the location of right click event for picking interogation
	intPair pickXY;
	
	// -1 if no object selected, else holds id of selected object
	int pickID;
	
	// holds the original colour and transparency value of objects
	vector<RGBA> RGBAHold;
	
	// flag, indicates whether the program is in scroll mode
	int scrollFlag;
	
	// holds the parameters for slideshow mode
	slideParam slideData;
	
	// flag that indicates whether slideshow mode is active
	int slideshowMode;
	
	// holds list of coi drawn when slideshow mode activated
	vector<intDraw> coiDrawHold;
	
	// holds list of objects that have key specified by user for slideshow mode
	vector<int> keySphereList; 
	
	// holds list of pairs for a key value if pair selection mode requested
	vector<intPair> keyPairList;
	
	// holds current position in keySphereList
	int keyListPos; 
	
	public slots:
	
	/**
	   slot to set colour of selected object
	*/
	void setColSlot(int colID);
	
	/**
	   slot to control changing the colour of selected objects
	*/
	void colChangeSlot();
	
	/**
	   slot to handle a change in the translider value
	*/
	void transSliderChanged(int val);
	
	/**
	   slot to control changing the transparency of selected objects
	*/
	void transChangeSlot();
	
	/**
	   slot that prints intersections for currently selected object to standard output
	   
	   @param id id of currently selected object
	*/
	void printIntersections(int id);
	
	/**
	   slot that updates coiDraw vector depending on currently selected object
	*/
	void updateDrawList();
	
	/**
	   slot that prints intersections for all objects to standard output
	*/
	void printAllIntersections();
	
	/**
	   slot that updates coiDraw vector for all objects
	   if all objects are drawn then clear coiDraw vector else add missing coi to coiDraw
	*/
	void drawAllIntersections();
	
	/**
	   Slot that creates the dialog that gets user parameters for slideshow mode
	*/
	void createSlideDialog();
	
	/**
	   Slot that sets up program for slideshow mode
	*/
	void setSlideshow();
	
	/**
	   Slot that handles  value changes in the key combo box
	*/
	void keyChange(QString key);
	
	/**
	   Slot that handles value changes in the transition combo box
	*/
	void transChange(QString trans);
	
	/**
	   Slot that handles value changes in the group combo box
	*/
	void groupChange(QString group);
	
	
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
	
	/**
	   interaction handling for when mouse button released, used for picking
	   
	   @param event information about the mouse button released
	*/
	void mouseReleaseEvent(QMouseEvent *event);
	
}; // end of class

#endif // ends our header guard
