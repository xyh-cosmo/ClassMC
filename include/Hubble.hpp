#ifndef __HUBBLE__
#define __HUBBLE__

#include <iostream>
#include <string>
#include <imcmc/imcmc.hpp>

#include "ParamList.hpp"
#include "Misc.hpp"

//  Hubble parameters measured using differential ages of passive evloving galaxies
/*  data format
    z   Hz  dHz
 */

struct Data_Hubble{

    DataInfo data_info;

    int     size;
    double  hubble_z_max;
    double  *z;
    double  *Hz;
    double  *dHz;
    bool    has_data;

    Data_Hubble();
    ~Data_Hubble();

	void Init( std::string& paramfile );
    void ReadData( std::string& paramfile );
    void MaxRedshift();
};

double Likelihood_Hubble( imcmc::imcmc_double&  param,
                          double&               lndet,
                          double&               chisq,
                          void*                 model,
                          void*                 data,
                          imcmc::istate&        state );

#endif  //  __HUBBLE__
