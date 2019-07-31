#include "CosmoTheory.hpp"
#include "SNE.hpp"

using namespace std;
using namespace imcmc;

double Likelihood_SNe_WFIRST(imcmc_double&  param,
                             double&        lndet,
                             double&        chisq,
                             void*          model,
                             void*          data,
                             istate&        state ){

    state.this_like_is_ok = true;
	lndet = chisq = 0;

    CosmoTheory	*cosmo	= static_cast<CosmoTheory*>(model);
	SNe_WFIRST 	*sne    = static_cast<SNe_WFIRST*>(data);

	double MB = param["MB"];
	double dmu;

	for( int i=0; i<sne->sne_num; ++i ){
		double dl = cosmo->engine->get_Dl(sne->z(i));
		dmu = 5.0*log10(dl)+ 25.0 - sne->mu(i) + (MB+19.3);
		chisq += pow( dmu/sne->dmu(i), 2 );
	}

	return -lndet - 0.5*chisq;
}
