/*    drivenMain.cpp

      Ben Postlethwaite
      Physics 410 Assignment 4

      *** Main ***
      */
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "drivenHeaders.h"

using std::endl; using std::ofstream;
using std::cout; using std::cerr;


int main(int numargs,char* args[])
{
  
  // ARGUEMENT CHECKING
  if (numargs != 4){
    cout << "You Need to enter 3 arguments following driven\n";
    cout << "./driven (Initial theta) (Initial Velocity) (Initial A)" << endl;
    return 0;
    }
  
  //  declarations of variables, calculating parameters
  double *y, *dydt, *yout;
  double tmax, tscale,t, omega0;
  double initial_theta = atof(args[1]), initial_v = atof(args[2]);
  //double nu = atof(args[3]),  A = 0;; // For part b
  double A = atof(args[3]) ,      nu = 0.5;   // for part c
  double PI = 3.14159265358;     
  double h = 0.01, l = 1,         g = 1;
  double omega, param[2],         m = 1;
  int i, numSteps, nequations = 2,nperiods = 300;
  omega = 2/3;
  omega0 = sqrt(g/l);
  param[0] = nu*omega0/(m*g);     //param is a list of parameters
  param[1] = A/(m*g);
  param[2] = omega/omega0;  

  // Open File ready for writing
  ofstream os;                    // os is like cin but to file
  const char * fn1 = "phys05.dat";
  // Make sure we opened the file
  os.open(fn1);
  if( !os ) {                      // file couldn't be opened
    cerr << "Error: file could not be opened" << endl;
    exit(1);
    } 
   
  // Create some arrays in memory
  dydt = new double[nequations];
  y    = new double[nequations];
  yout = new double[nequations];

  //  setting initial values, tmax and nondimensionalization
  tscale = sqrt(l/g);
  tmax = (2 * nperiods * PI) / tscale;// NON dimensionalize tmax
  initial_v = initial_v / tscale;
  y[0]   = initial_theta    ;        
  y[1]   = initial_v;                
  t=0.0;                             
  
  // now we start solving the differential equations using the RK4 method
  while (t <= tmax){
    derivatives(t, y, dydt,param);         // initial derivatives
    RK4(y,dydt,nequations,t,h,yout,param,derivatives);
    for (i = 0; i < nequations; i++) {
      y[i] = yout[i];
    }
    
    t += h/tscale;                   // NON dimensionalize t increment
    fsend(os,t, y, h);                  // write to file
  }

  // Delete arrays and close files
  delete [] y; delete [] dydt; 
  delete [] yout;
  os.close();
  
  
  return 0;
}   //  End of main function


