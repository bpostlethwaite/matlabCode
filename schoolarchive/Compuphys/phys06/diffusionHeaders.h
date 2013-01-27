/* diffusionHeaders.h
   Ben Postlethwaite

   ***  Headers  ***
*/

#ifndef diffusion_H_
#define diffusion_H_

#include <iosfwd>
using std::ostream;

// function prototypes


void fsend(ostream &, double,const int, double * );
void Explicit(const int ,const double ,const double * ,double * );
double uInitial(const double );
void tridiag(double *, double *, double *, double *, int);
void initmesh(double *, double *, double *, int, double);


// Const Variables


#endif
