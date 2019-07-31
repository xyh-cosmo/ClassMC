#ifndef __SNE_JLA__
#define __SNE_JLA__

#include <string>
#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>

#include "ParamList.hpp"

/////////////////////////////////////////////////////
//  JLA data structure is actually defined in jla.h
/////////////////////////////////////////////////////

double Likelihood_SNe_JLA(  imcmc::imcmc_double&    param,
                            double&                 lndet,
                            double&                 chisq,
                            void*                   model,
                            void*                   data,
                            imcmc::istate&          state );

struct SNe_JLA_Mock{
    DataInfo    data_info;
    int         sne_num;
    bool        has_err;
    bool        use_full_cov;
	bool		invert_err;	// used to check randomness of w transition
    arma::vec   z;
    arma::vec   mb, dmb;
    arma::mat   covmat_inv;
    SNe_JLA_Mock();
    ~SNe_JLA_Mock();
    void Init( std::string& dataset );
};

double Likelihood_SNe_JLA_Mock(	imcmc::imcmc_double&	param,
								double&					lndet,
								double&					chisq,
								void*					model,
								void*					data,
								imcmc::istate&			state );

#endif  //  __SNE__
