/*
 *	This is a unified interface to all the likelihood functions, because the Cosmological Model
 *	is shared by more than one likelihood functions, so put them together is more convenient.
 */

#include "ClassMC.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

double Likelihoods(	imcmc_double&   param,
                    double&         lndet,
                    double&         chisq,
                    void*           model,
                    void*           data,
                    istate&         state ){

    double lndet_temp = 0.0;
    double chisq_temp = 0.0;

//    Print_ImcmcDouble(param);
//    exit(0);

    CosmoTheory	*cosmo 	= static_cast<CosmoTheory*>(model);
    DataList	*dlist 	= static_cast<DataList*>(data);

//  update cosmology
    bool status = cosmo->update(param);
//    bool status = cosmo->engine->updateParValues(param);


    if( status ){
        cosmo->update_derived_params(param);  // derived parameters are updated only if the cosmo->update(param) successed !!
        state.this_like_is_ok = true;
        state.store_mesg(" ===== class successfully updated cosmology ===== ");
    }
    else{	//  if CLASS failed to update, then directly return
        state.this_like_is_ok = false;
        state.store_mesg(" ===== failed to update cosmology ===== ");
//        chisq = _IMCMC_CHISQ_MAX_;
        chisq = 1e100;
        return -0.5*chisq;
    }

    chisq = lndet = 0.0;

    if( dlist->use_Age ){
		Likelihood_Age( param,
                        lndet_temp,
                        chisq_temp,
                        cosmo,
                        dlist->data_age,
                        state );
		lndet += lndet_temp;
		chisq += chisq_temp;
    }

    if( dlist->use_BAO ){
     Likelihood_BAO( param,
                     lndet_temp,
                     chisq_temp,
                     cosmo,
                     dlist->data_bao,
                     state );
     lndet += lndet_temp;
     chisq += chisq_temp;
    }

    if( dlist->use_HST ){
		Likelihood_HST( param,
                        lndet_temp,
                        chisq_temp,
                        cosmo,
                        dlist->data_HST,
                        state );
		lndet += lndet_temp;
		chisq += chisq_temp;
    }

    if( dlist->use_Hz ){
		Likelihood_Hubble( 	param,
                            lndet_temp,
                            chisq_temp,
                            cosmo,
                            dlist->data_Hz,
                            state );
		lndet += lndet_temp;
		chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
cout << "Chi2_hubble = " << chisq_temp << endl;
#endif
    }

    if( dlist->use_RSD ){
        Likelihood_RSD( param,
                        lndet_temp,
                        chisq_temp,
                        cosmo,
                        dlist->data_rsd,
                        state );
        lndet += lndet_temp;
        chisq += chisq_temp;
    }

    if( dlist->use_SN ){
    	if( dlist->use_UNION ){
			Likelihood_SNe_UNION( 	param,
                            	lndet_temp,
                            	chisq_temp,
                            	cosmo,
                            	dlist->data_sne_union,
                            	state );
			lndet += lndet_temp;
			chisq += chisq_temp;
    	}
    	else if( dlist->use_SNLS ){
			Likelihood_SNe_SNLS(	param,
	                            lndet_temp,
	                            chisq_temp,
	                            cosmo,
	                            dlist->data_sne_snls,
	                            state );
			lndet += lndet_temp;
			chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
cout << "Chi2_snls = " << chisq_temp << endl;
#endif
    	}
        else if( dlist->use_SNLS_Mock ){
            Likelihood_SNe_SNLS_Mock(    param,
                                    lndet_temp,
                                    chisq_temp,
                                    cosmo,
                                    dlist->data_sne_snls_mock,
                                    state );
            lndet += lndet_temp;
            chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
cout << "Chi2_snls_mock = " << chisq_temp << endl;
#endif
        }
    	else if( dlist->use_JLA ){
//    		cout << "try calling JLA likelihood ...\n";
			Likelihood_SNe_JLA( param,
                            lndet_temp,
                            chisq_temp,
                            cosmo,
                            dlist->data_sne_jla,
                            state );
			lndet += lndet_temp;
			chisq += chisq_temp;

//			cout << "got JLA-chi2: " << chisq_temp << endl;
    	}
    	else if( dlist->use_JLA_Mock ){
    		Likelihood_SNe_JLA_Mock(param,
		                            lndet_temp,
		                            chisq_temp,
		                            cosmo,
		                            dlist->data_sne_jla_mock,
		                            state );
    		lndet += lndet_temp;
    		chisq += chisq_temp;
    	}


		if( dlist->use_WFIRST ){
			Likelihood_SNe_WFIRST(	param,
									lndet_temp,
									chisq_temp,
									cosmo,
									dlist->data_sne_wfirst,
									state );
			lndet += lndet_temp;
			chisq += chisq_temp;
		}

    }

    if( dlist->use_CMB ){
    /*
        Planck needs Cl,
        WMAP needs Cl*l*(l+1)/2/pi
    */
    	if( dlist->use_PLK ){
			Likelihood_PLK( param,
                            lndet_temp,
                            chisq_temp,
                            cosmo,
                            dlist->data_plk2015,
                            state );
			lndet += lndet_temp;
			chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
cout << "Chi2_plk15 = " << chisq_temp << endl;
#endif
    	}
#if defined (_USE_WMAP7_)
        else if( dlist->use_WMAP7 ){
            Likelihood_WMAP7( param,
                              lndet_temp,
                              chisq_temp,
                              cosmo,
                              dlist->data_wmap7,
                              state );
            lndet += lndet_temp;
            chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
cout << "Chi2_wmap7 = " << chisq_temp << endl;
#endif
        }
#elif defined (_USE_WMAP9_)
        else if( dlist->use_WMAP9 ){
            Likelihood_WMAP9( param,
                              lndet_temp,
                              chisq_temp,
                              cosmo,
                              dlist->data_wmap9,
                              state );
            lndet += lndet_temp;
            chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
cout << "Chi2_wmap9 = " << chisq_temp << endl;
#endif
        }
#endif
    }

    // cout << "chisq_tot = " << chisq << endl;

    //  if there are any priors
    if( dlist->use_DDE_CPZ_prior ){
        Prior_DDE_CPZ(  param,
                        lndet_temp,
                        chisq_temp,
                        cosmo,
                        dlist->prior_dde_cpz,
                        state );
        lndet += lndet_temp;
        chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
cout << "Chi2_dde = " << chisq_temp << endl;
#endif
    }

    if( dlist->use_CMB_DistPrior ){
    //  consitency with full CMB data should have been checked during the initialization
        Prior_CMB_Dist( param,
                        lndet_temp,
                        chisq_temp,
                        cosmo,
                        dlist->prior_cmb_dist,
                        state );
        lndet += lndet_temp;
        chisq += chisq_temp;
#ifdef _DEBUG_CHI2_
        cout << "Chi2_CMB_dist = " << chisq_temp << endl;
#endif
    }

#ifdef _DEBUG_CHI2_
cout << "******* Chi2_tot = " << chisq << endl;
#endif


    return lndet - 0.5*chisq;
}
