#include "CosmoTheory.hpp"
#include "SNE.hpp"

using namespace std;
using namespace imcmc;

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

// IMPORTANT: the absolute magnitude MB should be fit together with other cosmological parameters
double Likelihood_SNe_SNLS_Mock(imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state ){

    state.this_like_is_ok = true;

    lndet = chisq = 0;

    CosmoTheory     *cosmo  = static_cast<CosmoTheory*>(model);
    SNe_SNLS_Mock   *sne    = static_cast<SNe_SNLS_Mock*>(data);

    double MB = param["MB"];
    double dl_i, mu_i, mb_i;
    arma::vec dmb = arma::zeros(sne->sne_num,1);

	//never forget to fit MB, if it is not marginalized!!!
    for(int i=0; i<sne->sne_num; ++i){
        dl_i = cosmo->engine->get_Dl(sne->z(i));
        mu_i = 5.*log10(dl_i) + 25;
        mb_i = mu_i + MB;   // convert distance modulus to apparent magnitude
        dmb(i) = mb_i - (sne->mb(i) - 19.3); // sne->mb is actually distance modulus ...
    }

    chisq = arma::as_scalar( dmb.t() * sne->covmat_inv * dmb );

    // MPI_DEBUG_STOP;

    return -lndet - 0.5*chisq;
}
