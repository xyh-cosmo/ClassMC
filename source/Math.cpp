#include <iostream>
#include "math_libs.hpp"

// Romberg integrator:
// note this one is different from the one used inside CLASS
double Romberg_Integrator( 	double (*f)(double, void *),
                            double a, double b,
                            void *param,
                            double eps ){
    int m,n,i,j;
    double y[_ROMBERG_ITER_MAX_];
    double h,ep, p, x, s, q=-9999;

    h = b-a;
    y[0] = h*( (*f)(a,param) + (*f)(b,param) )/2.0;
    m = 1;
    n = 1;
    ep = eps + 1.0;

    while( (ep >= eps) && (m < _ROMBERG_ITER_MAX_) ){
        p = 0.0;
        for( i=0; i<=n-1; ++i ){
            x = a + (i+0.5)*h;
            p = p + (*f)(x,param);
        }

        p = (y[0] + h*p)/2.0;
        s = 1.0;

        for( j=1; j<=m; ++j ){
            s = 4.0*s;
            q = (s*p - y[j-1])/(s-1.0);
            y[j-1] = p;
            p = q;
        }

        ep = fabs(q-y[m-1]);
        m = m + 1;
        y[m-1] = q;
        n = n + n;
        h = h/2.0;
    }

    //	one possible problem is that even when the iteration number has exceeded _ROMBERG_ITER_MAX_,
    //	the precision requirement still has not been satisified
    if( (m >= _ROMBERG_ITER_MAX_) && (ep >= eps)){
        std::string err = "\n";
        err += "*** double Romberg_Integrator( double (*f)(double, void *), double a, double b, void *param, double eps ) ==>\n";
        err += "*** the iteration has exceeded the maximum number, but the presion still has not reached!";
    }

    return q;
}
