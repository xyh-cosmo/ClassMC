#ifndef __PLANCK2015__
#define __PLANCK2015__

#include <string>
#include <map>
#include <imcmc/imcmc.hpp>

#include "CosmoTheory.hpp"
#include "ParamList.hpp"

#include "clik.h"
#include "clik_lensing.h"

//  ==================================================
//  these two are placed here just for compilation
// struct CosmoTheory;
// struct Data_Planck2015;
//  ==================================================

using namespace imcmc;
using namespace imcmc::parser;

struct Data_Planck2015{

    DataInfo data_info;

    bool    has_highl;  // determined from clikfile_highl
    bool    has_lowP;   // determined from clikfile_lowP
    bool    has_lens;   // determined from clikfile_lens

//  ====    common variables    ====
    int             i,cli;
    error           *_err,**err;
    parname         clnames[6];

    int has_cl_highl[6], lmax_highl[6];
    int has_cl_lowP[6], lmax_lowP[6];
    int lmax_lens[7];  // lensing is differen from the above two

    int nextra_highl, ndim_highl;
    int nextra_lowP, ndim_lowP;
    int nextra_lens, ndim_lens;

    double loglike_highl;
    double loglike_lowP;
    double loglike_lens;
	double loglike_A_Planck;
    double loglike_total;

//  extra names
    parname *names_highl;
    parname *names_lowP;
    parname *names_lens;
    std::vector<std::string> extra_names;  // in case there are common variables in names_highl, names_lowP and names_lens (this is actually empty)

    double *cl_and_pars_highl;
    double *cl_and_pars_lowP;
    double *cl_and_pars_lens;

    clik_object         *highl;
    clik_object         *lowP;
    clik_lensing_object *lens;

    int nall_highl;
    int nall_lowP;
    int nall_lens;

    bool Init( std::string& dataset );
    void update_extra_params( imcmc_double& full_param );

    void get_cls( CosmoTheory* theory );
    void get_cls_highl( CosmoTheory* theory );
    void get_cls_lowP( CosmoTheory* theory );
    void get_cls_lens( CosmoTheory* theory );

    bool    has_A_Planck;   //  if ture, add chisq_{A_{Planck}}
    double  A_Planck_mean, A_Planck_std;    //  default value: 1.0 , 0.0025

    void compute_chisq();

    Data_Planck2015();
    ~Data_Planck2015();

    void get_Cl_tasks( bool& has_tCl, bool& has_pCl, bool& has_lCl );
};


double Likelihood_PLK(  imcmc_double&   full_params,
                        double&         lndet,
                        double&         chisq,
                        void*           model,
                        void*           data,
                        istate&         state );

#endif  //  __PLANCK2015__
