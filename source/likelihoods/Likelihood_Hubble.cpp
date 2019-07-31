#include "CosmoTheory.hpp"
#include "Hubble.hpp"

using namespace std;
using namespace imcmc;

double Likelihood_Hubble(   imcmc_double&    param,
                            double&          lndet,
                            double&          chisq,
                            void*            model,
                            void*            data,
                            istate&          state ){

    state.this_like_is_ok = true;

    lndet = chisq = 0;

    CosmoTheory *cosmo = static_cast<CosmoTheory*>(model);
    Data_Hubble *Hz = static_cast<Data_Hubble*>(data);

    for(int i=0; i<Hz->size; ++i){
        double zi   = Hz->z[i];
        double Hzi  = cosmo->engine->get_Hz_km_s_Mpc(zi);
        double dHz  = (Hzi - Hz->Hz[i]) * (Hz->z[i] <= Hz->hubble_z_max);
        chisq += pow(dHz/Hz->dHz[i], 2);
    }

    return -lndet - 0.5*chisq;
}
