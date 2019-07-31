//  ========================================================================
//  Author: Youhua Xu
//  Date:   Jan-7, 2017
//
//  Supernovae data structures and nuisance parameters
//
//  JLA is defined in another independent source file
//  ========================================================================

#ifndef __SNE_UNION__
#define __SNE_UNION__

#include <string>
#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>

#include "ParamList.hpp"

struct SNe_UNION{

    DataInfo    data_info;

    int         sne_num;    //  580
    arma::vec   z;          //  redshifts of SNe
    arma::vec   mu, dmu;    //  distance modulus and uncertainties
    arma::vec   P;          //  probability of low mass
    arma::mat   icov;       //  covariance matrix, two are provided, nosys / sys, has aready been inverted
    bool        systematic; //  true = sys, false = nosys

    SNe_UNION(){};  // do nothing
    ~SNe_UNION() {}; // do nothing

    double      max_z();
    void        Init( std::string& dataset );
};



double Likelihood_SNe_UNION(imcmc::imcmc_double&    param,
                            double&                 lndet,
                            double&                 chisq,
                            void*                   model,
                            void*                   data,
                            imcmc::istate&          state );



#endif  //  __SNE_UNION__
