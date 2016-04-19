/**
 * @author Cameron Kumar
 * @version 1.0 5/3/16
 * @date 5 Mar 2016
 * @brief class definition for 3d elliptic non-Euclidean visualisation tool.
 *
 * Class declaration for the visualisation part of the non-Euclidean 3D
 * visualisation software. An instance of this will appear within a QT GUI. The  
 * 3D graphics defined here is handled by OpenGL. Input is specified from the 
 * location this class is instanced.  
*/

#include <QGLWidget> 

using namespace std;

#ifndef __VIS_H__
#define __VIS_H__

/**
 * struct containing 3 double values for use as coordinate points.
 */  
struct Point {
	
	double x; /**< x coordinate of point. */
	double y; /**< y coordinate of point. */
	double z; /**< z coordinate of point. */
	
};

/**
 * structure to hold quads of floats for use as RGBA values
 */
struct RGBA {

	float R; /**< Red value of colour, between 0.0 and 1.0. */
	float G; /**< Blue value of colour, between 0.0 and 1.0. */
	float B; /**< Green value of colour, between 0.0 and 1.0. */
	float A; /**< Alpha value of colour, defines transparency, between 0.0 and 1.0. */

};

/**
 * structure to store generation and key data of objects.
 */
struct genKey {

	int gen; /**< positive integer for object's generation. */
	string key; /**< object's key value. */
	
};

/** a structure that contains an id and distance from plane value.
 *  used for sorting translucent objects.
 */
struct idDist {

	int id; /**< id value of object for which distance has been computed. */
	double dist; /**< distance from plane parallel to camera through origin. */

};

/**
 * structure that contains an id, overlap distance, direction vector and length.
 * used for sorting storing intersection details between two objects.
 */
struct idOverVecLen {
	
	int id; /**< id value of object that intersects object we are working with. */
	double over; /**< overlap, sum of their radii - distance between them. */
	Point vec; /**< direction vector from object we're working with's centre to id's centre. */
	double len; /**< distance between centre of both objects. */	
	
};

/**
 * structure to store data about a particular intersection, including circle of intersection details,
 * vector between the two spheres, and the ids of the two spheres involved. Used for drawing intersections.
 */
struct intDraw {

	Point cen; /**< centre of circle of intersection. */
	Point vec; /**< orthogonal vector to circle of intersection. */
	double rad; /**< radius of circle of intersection. */
	int id1; /**< id of one object involved in intersection. */
	int id2; /**< id of other object involved in intersection. */

};

/**
 * stores 2 integers. used when pair cannot be used ie. for non-compile-time constants
 */
struct intPair {
	
	int x; /**< first integer of pair. */
	int y; /**< second integer of pair. */
	
};

/**
 * stores data about the parameters of slideshow mode specified by the user
 */
struct slideParam {

	string key;  /**< key value of slideshow. */
	int trans;  /**< transition type, 0 if keystroke, 1 if time. */
	int group; /**< group type, 0 if individual, 1 if pair. */

};

/**
 * Class defining a Qt Widget that is an interactive visualisation to manipulate 
 * elliptic non-Euclidean geometry. Handles mouse and keyboard input with the 
 * use of protected functions inherited from QGLWidget. Data is passed to the 
 * visualisation from the program that instances it. 
 *
 * Information needs to be specified in text files with the each line corresponding 
 * to the details of one object in the following format,
 * 	x_centre_coordinate y_centre_coordinate z_centre_coordinate radius generation key
 * in which the centre coordinates and radius are real numbers, generation is an integer
 * and the key value is a string. 
 *
 * Camera Controls : left click and drag rotates camera, right click and drag moves camera,
 * scroll wheel or "-" and "+" keys zoom camera.
 * 
 * Scrolling the scroll wheel while an object is selected enters "selection mode". You may
 * use this mode to scroll through all spheres on the screen and press the right mouse 
 * button to select an object. The escape key will exit scroll mode without selecting an 
 * object
 *
 * Right mouse button will spawn a context menu, where information about intersections
 * can be output to the console or the intersections themselves can be drawn on 
 * screen for individual or all objects. Object's colour or transparency properties
 * can also be changed from this menu. Slideshow mode can be entered from here also,
 * where users can display objects of a specified key value one at a time or in pairs.
 * Slideshow mode may either transition to the next slide either when the space bar is 
 * pressed or after 10 seconds, depending on which option the user requests. Hitting the
 * escape key while in keypress mode will exit slideshow mode. There is also a "help" 
 * feature that will print this information to the console when selected.
 */
class vis: public QGLWidget {

	// Q_OBJECT macro
	Q_OBJECT

	public:
	
	/**
	 * class constructor
	 *  
	 * @param parent parent class, inherit from there 
	 */
	vis(QWidget *parent);
		
	/**
	 * creates the points of a sphere object, centre at origin, radius 1. this 
	 * template sphere can be copied and used to represent all spheres in our
	 * visualisation.
	 *  
	 * @param nPoints number of points per circle in the sphere
	 * @return vector of double triplets defining the point coords in 3d 
	 */
	vector<Point> makeSpherePoints(int nPoints);
	
	/**
	 * create points of a cube about the origin with side length 1.
	 * cube points added front to back, top to bottom, left to right.
	 * this cube will be transformed and drawn about objects when they are
	 * selected by the user.
	 *  
	 * @return vector containing points of cube
	 */
	vector<Point> makeCubePoints();
	
	/**
	 * draws a sphere with specified centre and radius.
	 *  
	 * @param centre coordinate point of centre location
	 * @param radius radius of sphere 
	 */
	void drawSphere(Point centre, double radius);
	
	/**
	 * writes centre, radius, key and generation information from file to 
	 * vectors. File specified as a command line argument. Sphere centre coordinates
	 * stored in objCentre vector, radius data in objRadius, and key and generation
	 * data in objGenKey.
	 * 
	 * @param objData file containing object raw data 
	 * @return returns 1 if there is a major error, else 0
	 */
	int setData(char* objData);
	
	/**
	 * populate the colours vector with a selection of RGBA values and creates colour
	 * icons from png files to be used in colour selection.
	 */
	void initColours();
	
	/**
	 * determines whether a string is already present within a list or not.
	 * 
	 * @param val value to be searched for
	 * @param list list to be searched in
	 * @return 1 if true, else 0
	 */
	int valPresent(string val, vector<string> list);
	
	/**
	 * finds all unique key values for all objects in the visualisation.
	 *  
	 * @return vector containing all unique key values
	 */
	vector<string> getKeyList();
	
	/**
	 * recursive function to be used for merge sorting of a vector.
	 *  
	 * @param vec vector to sort (structure of ints and doubles)
	 * @return sorted vector (structure of ints and doubles)
	 */
	vector<idDist> mergeSort(vector<idDist>);
	
	
	/**
	 * orders the spheres based on opacity and which is closest to the camera.
	 * this is done to order the spheres so they may be drawn from furthest to
	 * nearest to provide the see-through effect when rendered.
	 *  
	 * @return returns a vector of integers representing the order of object rendering.
	 */
	vector<int> sphereOrder();
	
	/**
	 * changes the RGB colour of object with index id to the RGB value provided.
	 *  
	 * @param id identifier of object whose colour we want to change
	 * @param rgb RGB colour we want the object to have
	 */
	void changeColour(int id, RGBA rgb);
	
	/**
	 * changes the transparency value of object to alpha value provided.
	 *  
	 * @param id identifier of object whose transparency we want to change
	 * @param alpha alpha value we want to change object to
	 */
	void changeTransparency(int id, double alpha);
	
	/**
	 * calculates the centre, radius, and orthogonal vector for an intersection.
	 * stores details in global vector coi. this information is later used to
	 * draw or print intersections.
	 * 
	 * @param id index of object to calculate intersections for
	 * @param inter list of intersecting objects
	 */
	void calculateIntersection(int id, vector<idOverVecLen> inter); 
	
	/**
	 * identifies which objects intersect specified object by comparing the
	 * sum of their radii with the distance between their centres.
	 *  
	 * @param id identifier of object we will calculate intersections for
	 * @return number of objects this object intersects with
	 */
	int intersectsWith(int id);
	
	/**
	 * draws a circle given a specified centre and radius. used to draw 
	 * intersections between spheres.
   	 *
   	 * @param cen centre of circle to draw
   	 * @param rad radius of circle to draw
	 */
	void drawCircle(intDraw circ);
	
	/** 
	 * draws the intersections saved in the global variable coiDraw.
	 */
	void drawIntersections();
	
	/**
	 * renders for picking using the colour hack. creates unique colour for 
	 * each object then determines which object the user has picked depending
	 * on the colour of the pixel under the mouse.
	 * 
	 * @return integer id of picked object
	 */
	int getPicked();
	
	/**
	 * draws cube with side lengths 1 about the origin from points held in
	 * cubePoints variable.
	 *  
	 * @param p array of int containing each point on cube face
	 */
	void drawCubeFace(int p[4]);
	
	/**
	 * renders a wireframe cube around currently selected picked object
	 *  
	 * @param id of sphere that cube should be drawn around
	 */
	void selectionCube(int id);
	
	/**
	 *  creates the context menu when picking occurs. available options depend
	 *  upon whether an object or the background is selected.
	 */
	void createContextMenu();
	
	/**
	 * scroll mode command to move selection to next object. either moves one
	 * index value forwards if wheel scrolled up or one id value backwards if
	 * wheel scrolled down. changes transparency values accordingly.
	 *  
	 * @param delta the rotation amount of the scroll wheel
	 */
	void scroll(float delta);
	
	/**
	 * function to print specific intersection between two objects and add interstion
	 * to draw list. used in slideshow pair mode.
	 *  
	 * @param a index of first sphere
	 * @param b index of second sphere
	 */
	void handleIntersection(int a, int b);
	
	/** 
	 *  resets all variables back to their original state when slideshow mode
	 *  is finished and resets flags.
	 */
	void stopSlideshow();
	
	/**
	 * function that makes to program wait for specified number of seconds.
	 *  
	 * @param s number of seconds to wait
	 */
	void waitFunc(int s);
	
	/**
	 * sets up objects so the next slide may be rendered in slideshow mode.
	 * sets transparencies and draws intersections.
	 */
	void createSlide();
	
	/**
	 * moves specified value from current position in vector to the front of the
	 * specified vector. u
	 *
	 * @param val specified value to be swapped to fron
	 * @param vec vector for swap to take place in
	 * @return reordered vector
	 */
	vector<int> bringToFront(int, vector<int>);
	
	private:
	
	/**
	 * define number of points in one circle of a sphere
	 */
	int circPoints;
	
	/**
	 * holds the coordinate points for a standard sphere
	 */
	vector<Point> spherePoints; 
	
	/**
	 * holds the coordinate points for standard cube
	 */
	vector<Point> cubePoints;
	
	/**
	 * holds the data of all geometric objects centres
	 */
	vector<Point> objCentre;
	
	/**
	 * holds the radius data of geometric objects
	 */
	vector<double> objRadius;
	
	/**
	 * holds the RGB colour values of each object (default red), and transparency value A (default 1.0)
	 */
	vector<RGBA> objColour;
	
	/**
	 * holds the generation and key value for each object
	 */
	vector<genKey> objGenKey;
	
	/**
	 * holds the RGB colour options for objects
	 */
	vector<RGBA> colour;
	
	/**
	 * vector to store all unqiue key valeus
	 */
	vector<string> key;
	
	/**
	 * holds the QIcon types for each colour to be used in combobox
	 */
	vector<QIcon*> colourIcon;
	
	/**
	 * holds the starting position when the mouse is clicked
	 */
	QPoint startPos;
	
	/** 
	 * holds the current zoom scale of the camera, used to scale translation
	 */
	double scaleFactor;
	
	/**
	 * holds current rotation about the yaw axis, used for calculating transforms
	 */
	double yRot;
	
	/**
	 * holds the current rotation about the pitch axis, used to prevent over rotation
	 */
	double pRot;
	
	/**
	 * holds the position that a particular objects cois begin at in the coi vector
	 */
	vector<int> coiBegin;
	
	/**
	 * holds information about circles of intersection
	 */
	vector<intDraw> coi;
	
	/**
	 * holds information about circles of intersection to be drawn 
	 */
	vector<intDraw> coiDraw; 
	
	/**
	 * flag to indicate whether a camera translation is currently taking place
	 */
	int trans;
	
	/**
	 * flag to indicate whether we wish to render for colour picking
	 */
	int colourPicking;
	
	/**
	 * stores the location of right click event for picking interogation
	 */
	intPair pickXY;
	
	/**
	 * -1 if no object selected, else holds id of selected object
	 */
	int pickID;
	
	/**
	 * holds the original colour and transparency value of objects
	 */
	vector<RGBA> RGBAHold;
	
	/**
	 * flag, indicates whether the program is in scroll mode
	 */
	int scrollFlag;
	
	/**
	 * holds the parameters for slideshow mode
	 */
	slideParam slideData;
	
	/**
	 * flag that indicates whether slideshow mode is active.
	 */
	int slideshowFlag;
	
	/**
	 * holds list of coi drawn when slideshow mode activated.
	 */
	vector<intDraw> coiDrawHold;
	
	/**
	 * holds list of objects that have key specified by user for slideshow mode.
	 */
	vector<int> keySphereList; 
	
	/**
	 * holds list of pairs for a key value if pair selection mode requested.
	 */
	vector<intPair> keyPairList;
	
	/**
	 * holds current position in keySphereList
	 */
	int keyListPos; 
	
	public slots:
	
	/**
	 * slot to set colour of selected object, called from colour dialog.
	 */
	void setColSlot(int colID);
	
	/**
	 * slot to control creation of the colour change dialog, spawned from
	 * context menu.
	 */
	void colChangeSlot();
	
	/**
	 * slot to handle a change in the translider value, called from transparency
	 * dialog.
	 */
	void transSliderChanged(int val);
	
	/**
	 * slot to control creation of the transparency change dialog, spawned from
	 * context menu.
	 */
	void transChangeSlot();
	
	/**
	 * slot that prints intersections for currently selected object to standard output.
	 * called from context menu.
	 *  
	 * @param id id of currently selected object
	 */
	void printIntersections(int id);
	
	/**
	 * slot that updates coiDraw vector depending on currently selected object.
	 * called from context menu.
	 */
	void updateDrawList();
	
	/**
	 * slot that prints intersections for all objects to standard output. called
	 * from context menu.
	 */
	void printAllIntersections();
	
	/**
	 * slot that updates coiDraw vector to include all objects. however if all 
	 * objects are currently drawn, clears coiDraw.
	 */
	void drawAllIntersections();
	
	/**
	 * slot that creates the dialog that gets user parameters for slideshow mode. 
	 * employs a grid layout within a box layout.
	 */
	void createSlideDialog();
	
	/**
	 * slot that sets up program for slideshow mode by setting variables and flags
	 * required, depending on user specified parameters. creates lists of which
	 * objects need to be highlighted for each slide and carries out entire 
	 * slideshow if time transition selected.
	 */
	void setSlideshow();
	
	/**
	 * slot that handles  value changes in the key combo box within the slideshow
	 * mode dialog.
	 *
	 * @param key newly selected string value of key
	 */
	void keyChange(QString key);
	
	/**
	 * slot that handles value changes in the transition combo box within the slideshow
	 * mode dialog.
	 *
	 * @param trans newly selected transition parameter
	 */
	void transChange(QString trans);
	
	/**
	 * slot that handles value changes in the group combo box within the slideshow
	 * mode dialog.
	 *
	 * @param group newly selected grouping parameter
	 */
	void groupChange(QString group);
	
	
	protected:
	
	/**
	 * initialises environment for OpenGL rendering when instance called. 
	 * determines which sphere geometry to use. initializes variables and flags.
	 * creates points for sphere and cube geometry. sets up camera and lighting.
	 */
	void initializeGL();
	
	/**
	 * changes size of viewport when widget resized
	 * 
	 * @param w the new width of the widget
	 * @param h the new height of the widget
	 */
	void resizeGL(int w, int h);
	
	/**
	 * draw a new frame. if standard rendering then a new frame is created based
	 * on the current variable settings of the program and swapped to front buffer.
	 * if picking render requested, objects rendered by getPicked function to 
	 * determine what object the user has selected.
	 */ 
	void paintGL();
	
	/**
	 * interaction handling for the mouses's scroll wheel, used for 
	 * camera zoom and scroll mode.
	 *  
	 * @param event information about the mouse button click
	 */
	void wheelEvent(QWheelEvent *event);
	
	/**
	 * interaction handling for when a button on the keyboard is pressed,  
	 * used for camera zooming and scroll mode, moving between slides in slideshow
	 * mode, and exiting scroll and slideshow mode.
	 *  
	 * @param event information about the key pressed
	 */
	void keyPressEvent(QKeyEvent *event);
	
	/**
	 * interaction handling for when a button on the mouse is pressed,  
	 * stores the start position of the mouse for translations and rotations
	 * 
	 * @param event information about the mouse button press
	 */
	void mousePressEvent(QMouseEvent *event);
	
	/**
	 * interaction handling for when mouse movements take place, used for
	 * camera translations and rotations.
	 *  
	 * @param event information about the mouse button press
	 */
	void mouseMoveEvent(QMouseEvent *event);
	
	/**
	 *  interaction handling for when mouse button released, used for picking
	 *  
	 *  @param event information about the mouse button released
	 */
	void mouseReleaseEvent(QMouseEvent *event);
	
}; // end of class

#endif // ends our header guard
