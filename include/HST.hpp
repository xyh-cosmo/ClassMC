#ifndef __HST__
#define __HST__

#include <string>
#include <imcmc/imcmc.hpp>
#include "ParamList.hpp"

struct Data_HST{

    DataInfo data_info;

    double H0, sigma_H0;
    void Init(std::string& paramfile);

    Data_HST();
    ~Data_HST();
};

double Likelihood_HST(  imcmc::imcmc_double&    param,
    double&                 lndet,
    double&                 chisq,
    void*                   model,
    void*                   data,
    imcmc::istate&          state );

#endif  //  __HST__
