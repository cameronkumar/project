/**
   Qt GUI and interactive OpenGL visualisation for non-Euclidean geometry   
   defined by the moebinv library in 3D space
   
   @author Cameron Kumar
   @version 1.0 5/3/16
*/

#include <QApplication> // this is a Qt application
#include "vis.h" // header defining the visualisation

/**
   main function for the whole library, gui defined here
   
   @param argc number of command line arguments
   @param argv command line arguments
   @return exit status of function, normal 0
*/
int main(int argc, char *argv[]) {
	
	QApplication app(argc, argv); // creating the application
	vis *visualisation = new vis(NULL); // visualisation instance
	
	// now to display the visualisation in app
	visualisation->show();
	visualisation->resize(1024,512);
	
	// run the application
	return app.exec();
	
}
