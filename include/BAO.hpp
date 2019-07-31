/**************************  BAO Likelihood needs the follwoing parameters:  ***********************
    double  H0;         //  Hubble constant
    double  Omega_b;    //  baryon density fraction
    double  Omega_c;    //  dark matter density fraction
    double  Omega_k;    //  curature

    also some derived parameters( Not to be sampled )
    double  Neff;       //  energy density ratio of relativistic nuetrino to photon, fixed to 3.04
    double  Omega_r;    //  radiation density, Omgear = Omegag*(1+0.2271*Neff);
    double  Omega_de;   //  cosmological constant / dark energy, Omegal = 1-Omega_m-Omega_k-Omega_r

    data_type:
        rs/DV   =   1
        A(z)    =   2
        DV/rs   =   3
        DV/Mpc  =   4
        DA/rs   =   5
        _bao_DR12_zhao_ = 6 #  (DA/rs, Hz*rs)
****************************************************************************************************/

#ifndef __BAO__
#define __BAO__

#include <vector>
#include <string>
#include <sstream>
#include <imcmc/imcmc.hpp>
#include "ParamList.hpp"
#include "Misc.hpp"

struct CosmoTheory;

/*
    This is the lowest level data structure for BAO.
*/

#ifndef __BAO_TYPE__
#define __BAO_TYPE__
    #define _bao_rs_over_DV_    1
    #define _bao_A_             2
    #define _bao_DV_over_rs_    3
    #define _bao_DV_            4
    #define _bao_DA_over_rs_    5
    #define _bao_DR12_zhao_     6
    #define _bao_DA_			7
#endif

struct Data_BAO_Base{

    int type;
    int size;       // normally this is just 1, but some BAO data provides more data,
                    // thus a covariance matrix is needed.
    double *z;      // effective redshifts
    double *val;    // 
    double *err;    // errors of the results
    double *icov;   // inverse covariance matrix

    int read_data(std::string& datafile);
    int read_icov(std::string& covfile);

    Data_BAO_Base();
    ~Data_BAO_Base();

    void print_type(){
        std::cout << "BAO data type is: ";
        switch(type){
            case _bao_rs_over_DV_:
                std::cout << "rs/Dv " << std::endl;
                break;
            case _bao_A_:
                std::cout << "A" << std::endl;
                break;
            case _bao_DV_over_rs_:
                std::cout << "Dv/rs" << std::endl;
                break;
            case _bao_DV_:
                std::cout << "Dv" << std::endl;
                break;
            case _bao_DA_over_rs_:
                std::cout << "DA/rs" << std::endl;
                break;
            case _bao_DR12_zhao_:
            	std::cout << "(DA/rs, Hz*rs)" << std::endl;
            	break;
            case _bao_DA_:
            	std::cout << "DA" << std::endl;
            	break;
            default:
                std::cout << "UNKNOWN!!!" << std::endl;
                break;
        }
    }
};

struct Data_BAO{

    DataInfo data_info;

    int bao_data_size;   // number of BAO data (a single data may have more than one BAO measurements)
    std::vector<Data_BAO_Base*> BAO_List;   //  store various BAO data

    Data_BAO();
    ~Data_BAO();

    bool use_Hu_fitting;    // use Wayne Hu fitting formula. default true
    void Init( std::string& bao_dataset );
};

// the above are implementd in source/likelihoods/Likelihood_BAO.cpp


void Compute_BAO_Vals(  CosmoTheory* cosmo,
                        double *bao_z,
                        double *bao_prediction,
                        int&    bao_size,
                        int&    bao_type,
                        bool&   use_Hu_fitting );

double Compute_BAO_Chi2(double *bao_data, 
                        double *bao_prediction, 
                        double *bao_err,
                        double *bao_icov,
                        int& n_bao);

double Likelihood_BAO(  imcmc::imcmc_double&    param,
                        double&                 lndet,
                        double&                 chisq,
                        void*                   model,
                        void*                   data,
                        imcmc::istate&          state );

#endif  //  __BAO__
