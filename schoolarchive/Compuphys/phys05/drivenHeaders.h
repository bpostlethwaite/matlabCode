/* drivenHeaders.h
   Ben Postlethwaite

   ***  Headers  ***
*/

#ifndef driven_H_
#define driven_H_

#include <iosfwd>
using std::ostream;

// function prototypes
void derivatives(double, double *, double *,const double param[]);
void getdata ( double&, double&);
void fsend(ostream &, double, double *, double);
void RK4(double *, double *, int, double, double,
	 double *,const double param[],
	 void (*)(double, double *, double *,const double param[]));
// Const Variables


#endif
