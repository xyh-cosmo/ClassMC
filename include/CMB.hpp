/*  calculate CMB related observables, depends on Cosmology calculator(s).
    
    Shift parameter: R = (1+z*) * DA(z�6�5) * sqrt(Omegam) * H0 / c
    acoustic scale: lA = (1+z*) * Pi * DA(z*) / rs(z*)

where z* is the decoupling redshift
*/

#ifndef __CMB_DISTANCE_PRIOR__
#define __CMB_DISTANCE_PRIOR__

#include <vector>
#include <string>
#include <armadillo>
#include <imcmc/imcmc.hpp>
#include "ParamList.hpp"

using namespace imcmc;

struct Data_CMB_Dist {

    DataInfo data_info;
    bool use_Hu_fitting;
    int format;                     // 0--WMAP; 1--Planck
    arma::rowvec distance_prior;    // {lA,R,zdec}
    arma::mat covmat_inv;           // inverse of the covariance matrix.
    void Init( std::string& CMB_dist_prior_dataset );

    Data_CMB_Dist();
    ~Data_CMB_Dist();

//  add Planck CMB distance prior
    arma::rowvec distance_prior_plk;    // {R,lA,Omegabh2}
    arma::mat covmat_inv_plk;           // this (inv-)covariance matrix is different from those from WMAP
    double std_R, std_lA, std_Obh2;     // 1-sigma errors for {R,lA,Omegabh2}
};

typedef Data_CMB_Dist CMB_Dist;

//  likelihood function prototype:
double Prior_CMB_Dist(  imcmc_double&   param,
                        double&         lndet,
                        double&         chisq,
                        void*           model,
                        void*           data,
                        istate&         state );

#endif  //  __CMB_DISTANCE_PRIOR__
