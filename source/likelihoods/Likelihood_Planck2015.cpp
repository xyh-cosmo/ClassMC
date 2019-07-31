//    ==========================================================================
//    data is Planck CMB power spectrum
//    Each time when full_params is updated, class will calculate the Cls again,
//    after that clik will compute the corresponding likelihood
//    ==========================================================================

#include "ClassMC.hpp"

double Likelihood_PLK(  imcmc_double&   full_params,
                        double&         lndet,
                        double&         chisq,
                        void*           model,
                        void*           data,
                        istate&         state ){

    lndet = chisq = 0.0;

    CosmoTheory*        theory = static_cast<CosmoTheory*>(model);
    Data_Planck2015*    cldata = static_cast<Data_Planck2015*>(data);

    cldata->get_cls(theory);
    cldata->update_extra_params(full_params);
    cldata->compute_chisq();
    chisq = -2.0*cldata->loglike_total;

    return cldata->loglike_total;    //    imcmc always needs lndet - 0.5*chi^2
}
