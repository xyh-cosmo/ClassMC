#ifndef _F_SIGMA8_HPP_
#define _F_SIGMA8_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <imcmc/imcmc.hpp>
#include "ParamList.hpp"
#include "Misc.hpp"

struct CosmoTheory;

struct Data_RSD{
    DataInfo data_info;
    
    int size;
    double *zeff;   // effective redshifts
    double *val;    // RSD measurements
    double *err;    // errors
    // double *icov;   // no covariance currently

    Data_RSD();
    ~Data_RSD();

    int Init(std::string& dataset);
    int Read_Data(std::string& datafile);
};


void Compute_RSD_Vals(  CosmoTheory* cosmo,
                        double* zeff,
                        double* rsd,
                        int&    size);

double Likelihood_RSD(  imcmc::imcmc_double&    param,
                        double&                 lndet,
                        double&                 chisq,
                        void*                   model,
                        void*                   data,
                        imcmc::istate&          state );

#endif