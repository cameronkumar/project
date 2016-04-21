/**
   Qt GUI and interactive OpenGL visualisation for non-Euclidean geometry   
   defined by the moebinv library in 3D space
   
   @author Cameron Kumar
   @version 1.0 5/3/16
*/

#include <QApplication> // this is a Qt application
#include "vis.h" // header defining the visualisation
#include <iostream>

using namespace std;

/**
 * main function. sets up the Qt context, checks whether a file has been specified 
 * as a command line argument and creates an instance of the openGL visualisation.
 *  
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @return exit status of function, normal 0
 */
int main(int argc, char *argv[]) {
	
	QApplication app(argc, argv); // creating the application
	
	// want to read sphere data from command line argument, handled here
	if(argc == 2) { // create vis if 1 file specified
	
		vis *visualisation = new vis(NULL); // create vis instance
		
		// load data from file into visualisation
		if((visualisation -> setData(argv[1])) == 1) // error case
			return 1;
		else {
			// display the visualisation in app
			visualisation -> show();
			visualisation -> resize(1024,512);
	
			// run the application
			return app.exec();
		}
	} 
	
	else {
		cout << "ERROR: no geometric data file specified!\n";
		return 1;
	}	
	
	return 0;
}
