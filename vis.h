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
	vector<double> makeSpherePoints(int nPoints);

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

}; // end of class

#endif // ends our header guard
