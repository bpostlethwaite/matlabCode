/*  drivenDefs.cpp
    Ben Postlethwaite
    Physics 410 

    *** Function Definitions  ***

 */

// Function Declarations

//#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "diffusionHeaders.h" // Include protypes

using std::cout; using std::cin; using std::ios;
using std::endl; using std::setw; using std::setprecision;


//**********************************************************

//    function to write out the final results
void fsend(ostream & os,double t,const int n_mesh, double * u)
{
  int i;
  os << setw(15) << setprecision(4) << t;
  for(i=0;i<n_mesh+1;i++){
    os << setw(15) << setprecision(8) << u[i];
  }
   os << endl;
} 


//**************************************************************

// Create initial conditions:
double uInitial(const double x)
{
  return (0.0*x);
  }

//**************************************************************
// Explicit PD routine
void Explicit(const int i,const double alpha ,const double * u ,double * u2)
{
  u2[i] = alpha * u[i-1] + (1 - 2*alpha) * u[i] + alpha * u[i+1];
}
// **************************************************************
void initmesh(double *l, double *c, double *r, int N, double alpha)
{
  int i;
  for(i=1; i< N; i++)
    {
      l[i] = -alpha;
      c[i] = 1 + 2*alpha;
      r[i] = -alpha;
    }
  l[0]   = 0.0;
  l[N-1]   = 0.0;
  r[0]   = 0.0;
  r[N-1]   = 0.0;
}

//*****************************************************************

void tridiag(double *l, double *c, double *r, double *u2, int N)
{
  int i;
  for (i = 1; i < N-1;i++)
      {
	c[i]  = c[i] - l[i]*r[i-1];
	u2[i] = (u2[i] - l[i]*u2[i-1])/c[i];
	r[i]  = r[i] / c[i];
      }
  for ( i= N-1; i > 1;i--)
    {
      u2[i-1] = -r[i-1]*u2[i] + u2[i-1];
    }
}
