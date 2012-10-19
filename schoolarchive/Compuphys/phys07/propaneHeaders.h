/* propaneHeaders.h
   Ben Postlethwaite

   ***  Headers  ***
*/

#ifndef propane_H_
#define propane_H_

#include <iosfwd>
using std::ostream;

// function prototypes


void fsend(ostream &, double,const int, double *, double * );
void Explicit(const double,const double,const int, 
	      const double *,double *,double *, double *,
	      double, double *, double);
void initialize(double *, double *,double *, double *,
		double *, const int, const double,
		double, double, double *,double *);
void tridiag(double *, double *, double *, double *, int);
void initmesh(double *, double *, double *, double *,
	      double *, double, int, double, double *);


// Const Variables


#endif
