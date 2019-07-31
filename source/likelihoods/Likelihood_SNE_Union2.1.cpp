#include "CosmoTheory.hpp"
#include "SNE.hpp"

using namespace std;
using namespace imcmc;

//  =======================================================================
//  Union2.1 data is taken from cosmomc, in which H0 has been marginalized.
//  You may check that change H0 will not change chi2 ...
//  =======================================================================

double Likelihood_SNe_UNION( imcmc_double&  param,
                             double&        lndet,
                             double&        chisq,
                             void*          model,
                             void*          data,
                             istate&        state ){

    state.this_like_is_ok = true;
	lndet = chisq = 0;

    CosmoTheory	*cosmo	= static_cast<CosmoTheory*>(model);
	SNe_UNION 	*sne    = static_cast<SNe_UNION*>(data);

	arma::vec dmu = arma::zeros(sne->sne_num);

	for( int i=0; i<sne->sne_num; ++i ){
		double dl = cosmo->engine->get_Dl(sne->z(i));
		dmu(i) = 5.0*log10(dl)+ 25.0 - sne->mu(i);
	}

	chisq = arma::as_scalar( dmu.t() * sne->icov * dmu );

	return -lndet - 0.5*chisq;
}
