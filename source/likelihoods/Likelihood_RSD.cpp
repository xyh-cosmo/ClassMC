
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "CosmoTheory.hpp"
#include "RSD.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

void Compute_RSD_Vals(  CosmoTheory* cosmo,
                        double* zeff,
                        double* rsd,
                        int&    size){
    for( int i=0; i<size; ++i ){
        double f, sigma8;
        f       = cosmo->engine->get_f(zeff[i]);
        sigma8  = cosmo->engine->get_sigma8(zeff[i]);
        rsd[i]  = f*sigma8;

// for debug
//         cout << "## z = " << std::setw(7) << zeff[i]
//              << "f = " << std::setw(10) << f
//              << std::setw(10) << " sigma8 = " << sigma8 << endl;

//        if( fabs(sigma8) < 1e-5 ){
//            double Om = cosmo->engine->get_Omega_m();
//            cout << "#bug ==> Om = " << Om << endl;
//            exit(0);
//        }
    }

    // exit(0);
}

double Likelihood_RSD(  imcmc::imcmc_double&    param,
                        double&                 lndet,
                        double&                 chisq,
                        void*                   model,
                        void*                   data,
                        imcmc::istate&          state ){

    lndet = chisq = 0.0;

// for debug
//    imcmc_double_iterator it;
//    it = param.begin();
//    while( it != param.end() ) {
//        cout << it->first << " = " << it->second << "\n";
//        ++it;
//    }


    CosmoTheory	*cosmo  = static_cast<CosmoTheory*>(model);
    Data_RSD    *rsd 	= static_cast<Data_RSD*>(data);

    double *rsd_prediction = new double[rsd->size];
    Compute_RSD_Vals(	cosmo,
                        rsd->zeff,
                        rsd_prediction,
                        rsd->size );
    for( int i=0; i<rsd->size; ++i ) {
//         cout << "z = " << setw(10) << rsd->zeff[i]
//              << "prediction: " << setw(10) << rsd_prediction[i]
//              << "obs: " << setw(10) << rsd->val[i] << endl;

//        if( fabs(rsd_prediction[i]) < 1e-5 ){
//            fail = true;
//        }
        
        chisq += pow((rsd_prediction[i]-rsd->val[i])/rsd->err[i],2);
    }
//    cout << "-------------------------\n";
//     exit(0);    


    delete[] rsd_prediction;

    return lndet - 0.5*chisq;
}
