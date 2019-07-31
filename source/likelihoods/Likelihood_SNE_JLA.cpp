#include "CosmoTheory.hpp"
#include "SNE.hpp"
#include "jla.h"

using namespace std;
using namespace imcmc;


//  ==================================================================================================
//	JLA likelihood function wrapper (v4)
//  NOTE: JLA handles intrinsic dispersions (as well as sigma_z and sigma_lens) differently from SNLS3
//  More details can be found in the head of supernovae_JLA.f90 of cosmomc
//  ==================================================================================================
double Likelihood_SNe_JLA(  imcmc_double&   param,
                            double&         lndet,
                            double&         chisq,
                            void*           model,
                            void*           data,
                            istate&         state ){

    state.this_like_is_ok = true;

	lndet = chisq = 0;

//  ============================================================================================
//  update model parameters
//  by default, pass a CosmologyWorkspace object
//  ============================================================================================
	CosmoTheory	*cosmo	= static_cast<CosmoTheory*>(model);
	SNe_JLA 	*sne    = static_cast<SNe_JLA*>(data);

//  ============================================================================================
//  update LC nuisance parameters
//  NOTE: JLA::UpdateCov() actually doing nothing, all real calculations are done inside JLA
//        code provided by "Marc Betoule"
//  ============================================================================================
	sne->UpdateNuisance(param);
	sne->UpdateCov();

	double nuisance[4] = { sne->alpha, sne->beta, sne->MB, sne->DeltaMB };

	double *z   = sne->getZ();
	double *mu	= new double[sne->size()];

	double H0 = cosmo->engine->get_H0();
	double log10_70_over_H0 = log10(70./H0);

	for( int i=0; i<sne->size(); ++i ){
		mu[i] = 5.0*log10(cosmo->engine->get_Dl(z[i])) + 25.0 - 5.0*log10_70_over_H0;
		// mu[i] = 5.0*log10(cosmo->engine->get_Dl(z[i])) + 25.0;
	}

	chisq = sne->computeLikelihood(mu, nuisance);

	delete[] mu;
	delete[] z;

	return -lndet - 0.5*chisq;
}


double Likelihood_SNe_JLA_Mock( imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state ){
	state.this_like_is_ok = true;
	lndet = chisq = 0.0;

	//	========================
	//	update model parameters
	CosmoTheory		*cosmo	= static_cast<CosmoTheory*>(model);
	SNe_JLA_Mock	*sne	= static_cast<SNe_JLA_Mock*>(data);

	double MB = param["MB"];	// absolute peak magnitude in B band.
	double dl_i,mu_i,mb_i;

	arma::vec dmb = arma::zeros(sne->sne_num,1);
	//	========================================================================
	//	mock data use mb(z) as observables, this is different from real JLA data
	for( int i=0; i<sne->sne_num; ++i ){
		dl_i = cosmo->engine->get_Dl(sne->z(i));
		mu_i = 5.0*log10(dl_i) + 25.0;
		mb_i = mu_i + MB;
		dmb(i) = mb_i - (sne->mb(i) - 19.3);	// sne->mb is actually distance modulus
	}

	chisq = arma::as_scalar( dmb.t() * sne->covmat_inv * dmb );


	//exit(0);
	return -lndet - 0.5*chisq;
}
