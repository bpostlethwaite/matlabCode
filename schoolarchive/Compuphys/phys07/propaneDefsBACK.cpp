/*  propaneDefs.cpp
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
#include "propaneHeaders.h" // Include protypes

using std::cout; using std::cin; using std::ios;
using std::endl; using std::setw; using std::setprecision;


//**********************************************************

//    function to write out the final results
void fsend(ostream & os,double t,const int nmesh, double * u)
{
  int i;
  os << setw(15) << setprecision(4) << t;
  for(i=0;i<nmesh+1;i++){
    os << setw(15) << setprecision(8) << u[i];
  }
   os << endl;
} 


//**************************************************************
// initialize matrix
void initialize(double *u, double *u2,double *R, double *Rhalf,
		double *S, const int nmesh,const double deltaR,
		double K, double delta_t)
{  
  for (int i = 1; i <= nmesh; i++) 
    {
    u[i]  = 0;
    u2[i] = 0;
    R[i]  = i*deltaR;
    S[i]  = 0;
    }
  R[0] = 0;
  S[0] = S[1] = S[2] = S[3] = 1*K*delta_t;
  for (int i = 1; i <= nmesh; i++){
    Rhalf[i] = 0.5*(R[i-1] + R[i]);
  }
  u[0] = u[nmesh] = 4;
}
//**************************************************************
// Explicit PD routine
void Explicit(const double alpha ,const double K,const int nmesh,
	      const double * u ,double * u2,double *R, double *Rhalf)
{
for (int i = 1; i <= nmesh; i++) {
  u2[i] = alpha*K*Rhalf[i+1]*u[i+1]/R[i] + (1-2*alpha*K)*u[i] + 
    alpha*K*Rhalf[i]*u[i-1]/R[i];
 }
 u2[0]     = alpha*K*u[1] + (1-alpha*K)*u[0];
 u2[nmesh] = 0;
}
// **************************************************************
void initmesh(double *l, double *c, double *r, double *R, double *Rhalf,
	      double K, int N, double alpha, double delta_t, double *S)
{
  int i;
  for(i=1; i< N; i++)
    {
      l[i] = -K*alpha*Rhalf[i]/R[i] + S[i]*K*delta_t;
      c[i] = 1 + 2 * K * alpha + S[i]*K*delta_t;
      r[i] = -K*alpha*Rhalf[i+1]/R[i] + S[i]*K*delta_t;
    }
  l[0]   = 0;
  l[N-1]   = 0.0;
  r[0]   = 0;
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
