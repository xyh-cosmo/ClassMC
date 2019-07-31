#ifndef __DATA_AGE__
#define __DATA_AGE__

#include <vector>
#include <string>
#include <imcmc/imcmc.hpp>
#include "ParamList.hpp"

struct Data_Age {

    DataInfo data_info;

    int 	age_num;
    bool	age_init;
    double	*z;
    double	*age, *age_err;

    double	gft;		//	Galaxy formation time
    bool	gft_marg;	//	used to control which chisq will be used, default = false

    void Init( std::string& age_dataset );
    void Read_Data( std::string& data_file );

    Data_Age();
    ~Data_Age();
};

double Likelihood_Age(  imcmc::imcmc_double&    param,
                        double&                 lndet,
                        double&                 chisq,
                        void*                   model,
                        void*                   data,
                        imcmc::istate&          state );

#endif	//	__DATA_AGE__
