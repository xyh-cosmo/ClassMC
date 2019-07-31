#include "CosmoTheory.hpp"
#include "SNE.hpp"
#include "jla.h"

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

double Likelihood_SNe_SNLS( imcmc_double&   param,
                            double&         lndet,
                            double&         chisq,
                            void*           model,
                            void*           data,
                            istate&         state ){

    state.this_like_is_ok = true;

	lndet = chisq = 0;

    CosmoTheory	*cosmo	= static_cast<CosmoTheory*>(model);
	SNe_SNLS 	*sne    = static_cast<SNe_SNLS*>(data);

// 	update Light-Curve nuisance parameters: {alpha, beta}
//  MB is needed if sne->marg_scriptm == false
	sne->UpdateNuisance(param);
	sne->UpdateCov();

	arma::vec mb   = arma::zeros(sne->sne_num);;

	double dl_hel;
	double H0 = cosmo->engine->get_H0();	  // unit: km/s/Mpc
    double H0_c = H0 / (1e-3*_c_);
	double dl_with_cH0;

    for(int i=0; i<sne->sne_num; ++i){
        dl_with_cH0 = cosmo->engine->get_Dl(sne->zcmb(i));
        dl_hel	= dl_with_cH0 * H0_c;        // remove c_H0 dependence
        dl_hel	= dl_hel * (1.0 + sne->zhel(i)) / (1.0 + sne->zcmb(i));
        mb(i)   = 5.*log10(dl_hel);
    }

    if( sne->marg_scriptm ){

        double A,B,C,D,E,F,G;
        double estimated_scriptm, wtval;

        arma::vec invvars = arma::zeros(sne->sne_num,1);

        for( int i=0; i<sne->sne_num; ++i ){
            invvars(i) = 1.0 / (sne->pre_vars(i)
                + sne->alpha*sne->alpha*sne->ds(i)*sne->ds(i)
                + sne->beta*sne->beta*sne->dc(i)*sne->dc(i)
                + 2*sne->alpha*sne->cov_m_s(i)
                - 2*sne->beta*sne->cov_m_c(i)
                - 2*sne->alpha*sne->beta*sne->cov_s_c(i) );
        }

        wtval = arma::sum(invvars);
        estimated_scriptm = arma::as_scalar( (mb-sne->mb).t()*invvars ) / wtval;

        arma::vec dmb = sne->mb - mb
                        + sne->alpha*(sne->s-1)
                        - sne->beta*sne->c
                        - estimated_scriptm*0.;

        invvars = sne->icov_tot*dmb;
        A = arma::as_scalar(invvars.t()*dmb);

        if( sne->twoscriptmfit ){
            B = arma::as_scalar( invvars.t()*sne->K1 );
            C = arma::as_scalar( invvars.t()*sne->K2 );

            invvars = sne->icov_tot*sne->K1;
            D = arma::as_scalar( invvars.t()*sne->K2 );
            E = arma::as_scalar( invvars.t()*sne->K1 );

            invvars = sne->icov_tot*sne->K2;
            F = arma::as_scalar( invvars.t()*sne->K2 );

            G = E*F-D*D;
            chisq = A + log(0.5*E/_PI_) + log(0.5*G/_PI_/E) - B*B*F/G - C*C*E/G + 2.0*B*C*D/G;
        }
        else{
            arma::vec unit_colvec = arma::ones(sne->sne_num,1);
            B = arma::as_scalar(dmb.t()*sne->icov_tot*unit_colvec);
            E = arma::as_scalar(unit_colvec.t()*sne->icov_tot*unit_colvec);
//			E = arma::accu(sne->icov_tot);
            chisq = A + log(0.5*E/_PI_) - B*B/E;
        }
    }
    else{

        if( sne->twoscriptmfit ){
            for(int i=0; i<sne->sne_num; ++i)
                mb(i) +=  (-sne->alpha*(sne->s(i)-1) + sne->beta*sne->c(i) )
                        + sne->scriptMB1 + ( sne->var3(i) > sne->scriptmcut )*(sne->scriptMB2-sne->scriptMB1);
        }
        else{
            for(int i=0; i<sne->sne_num; ++i)
                mb(i) +=  (-sne->alpha*(sne->s(i)-1) + sne->beta*sne->c(i) )
                        + sne->scriptMB;
        }

        arma::vec dmb = sne->mb - mb;
        chisq = arma::as_scalar( dmb.t() * sne->icov_tot * dmb );
    }

	return -lndet - 0.5*chisq;
}

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
