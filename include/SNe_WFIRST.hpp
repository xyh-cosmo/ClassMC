/* format of the mock WFIRST sample
	z   mu+err   std    mu

There are 4 cols in total, where the last column is the fiducial value, which can
be used to detect biases in parameter constraints.
*/

#ifndef __SNE_WFIRST__
#define __SNE_WFIRST__

#include <string>
#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>

#include "ParamList.hpp"

struct SNe_WFIRST{

    DataInfo    data_info;

	bool		has_err;	// if true, then random errors are included.
    bool        invert_err; // if true, then multiply the random errors by -1

    int         sne_num;    //  number of SNe Ia
    arma::vec   z;          //  redshifts of SNe
    arma::vec   mu, dmu;    //  distance modulus and uncertainties

    SNe_WFIRST();
    ~SNe_WFIRST();

    void        Init( std::string& dataset );
};



double Likelihood_SNe_WFIRST(imcmc::imcmc_double&    param,
                            double&                 lndet,
                            double&                 chisq,
                            void*                   model,
                            void*                   data,
                            imcmc::istate&          state );



#endif  //  __SNE_WFIRST__
