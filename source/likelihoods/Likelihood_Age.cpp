#include "ClassMC.hpp"

using namespace std;
using namespace imcmc;

double Likelihood_Age( imcmc_double& 	param,
                       double& 			lndet,
                       double& 			chisq,
                       void* 			model,
                       void* 			data,
                       istate&          state ) {

    state.this_like_is_ok = true;

    lndet = chisq = 0;

    cout << " Caution: This likelihood function has not been finished !!!\n";

//    CosmoTheory	*cosmo  = static_cast<CosmoTheory*>(model);
//    Data_Age	*age 	= static_cast<Data_Age*>(data);

//    double dage;
//    double *age_theory = new double[age->age_num];
//    Age_at_z( age->z, age_theory, age->age_num, cosmo );

//    if( age->gft_marg ){	// use gft marginalized chisq

//    	for( int i=0; i<age->age_num; ++i ){
//    		dage = age->age[i] - age_theory[i];
//    		chisq += gsl_sf_log_erfc(dage/age->age_err[i]/sqrt(2.));
//    	}

//    	chisq *= (-2.0);
//    }
//    else{	// use nomral chisq, including gft as a model parameter

//    	age->Params.UpdateParams(param);
//    	age->gft = age->Params["gft"];
//    	double gft = age->gft;

//    	for( int i=0; i<age->age_num; ++i ){
//    		dage = age_theory[i] - (age->age[i] + gft);
//    		chisq += pow(dage/age->age_err[i], 2);
//    	}
//    }

    return -lndet - 0.5*chisq;

}
