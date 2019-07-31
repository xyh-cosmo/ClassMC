#include "HST.hpp"

using namespace std;
using namespace imcmc;

double Likelihood_HST( imcmc_double&	param,
                    	double&         lndet,
                    	double&         chisq,
                    	void*           model,
                    	void*           data,
                        istate&         state ){

    lndet = chisq = 0;

    Data_HST* HST = static_cast<Data_HST*>(data);

    chisq = pow((param["H0"] - HST->H0)/HST->sigma_H0, 2);

    return -lndet - 0.5*chisq;
}
