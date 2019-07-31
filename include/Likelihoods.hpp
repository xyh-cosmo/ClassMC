//	========================================================================
//	Author: Youhua Xu
//	Date:	Jan-7, 2017
//
// 	The likelihood function defined here is actually the "sum" of all other
// 	likelihood functions.
//
//	Individual likelihood functions are defined inside their corresponding
//	data structure headers.  By doing so the structure of ClassMC becomes
//	more simple and clear.
//	========================================================================


#ifndef __LIKELIHOODS__
#define __LIKELIHOODS__

#include <imcmc/imcmc.hpp>
#include <imcmc/ensemble.hpp>

#include "CosmoTheory.hpp"
#include "DataList.hpp"

using namespace imcmc;

double Likelihoods( imcmc_double&   full_params,
                    double&         lndet,
                    double&         chisq,
                    void*           model,
                    void*           data,
                    istate&         state );

// add all used likelihoods and data to ensemble_workspace
//  but use this has some unsolved problem ...
void Add_Likelihoods(   ensemble_workspace&         ew, 
                        CosmoTheory&                cosmo, 
                        DataList&                   datalist,
                        std::vector<string>&        mcmc_params,
                        std::vector<string>&        derived_params );

#endif // __LIKELIHOODS__
