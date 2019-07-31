#ifndef __MATH_LIBS__
#define __MATH_LIBS__

//	some constants

#ifndef _pi_
    #define _pi_    	3.14159265358979323846264338328
#endif

#ifndef _10_log10_
	#define _10_log10_ 	23.025850929940461
#endif


#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#define _ROMBERG_ITER_MAX_ 100

double Romberg_Integrator( 	double (*f)(double, void *),
                            double a, double b,
                            void *param,
                            double eps );

#endif
