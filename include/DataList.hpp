//	========================================================================
//	Author: Youhua Xu
//	Date:	Jan-7, 2017
//
// 	Data structure that contains all data, CMB, SNe, etc ...
//
//  Priors are also considered as data ...
//
//  Note: the structures of organizing datasets and likelihood functions can
//  be further optimized, for example, define a new "datatype", which pairs
//  each likelihood function to a specific data, and together the two parts
//  are stored in a std::vector constainor, so that the "if-else" would be
//  unnecessay, and the running time will also be decreased.
//	========================================================================
#ifndef __DATALIST__
#define __DATALIST__

#include "Age.hpp"
#include "BAO.hpp"
#include "HST.hpp"
#include "Hubble.hpp"
#include "CMB.hpp"
#include "RSD.hpp"

#if defined(_USE_WMAP7_)
    #include "WMAP7.hpp"
#elif defined(_USE_WMAP9_)
    #include "WMAP9.hpp"
#endif
#include "Planck2015.hpp"

#include "SNE.hpp"
#include "jla.h"

#include "ParamList.hpp"
#include "Prior.hpp"

#include "Misc.hpp"

using namespace imcmc;
using namespace imcmc::parser;

struct DataList{

//	=============
//	Data Options:
//	=============

    ParamList *ParList;  // use to collect MCMC & derived parameters

	bool use_Age;
	bool use_BAO;
    bool use_HST;       // HST H0 measurements
    bool use_Hz;        // Hubble parameters derived from passive evolving galaxies
    bool use_RSD;

    bool use_CMB;
#if defined(_USE_WMAP7_)
    bool use_WMAP7;		// WMAP7
#elif defined(_USE_WMAP9_)
    bool use_WMAP9;     // WMAP9
#endif
	bool use_PLK;       // Planck CMB data
	
	bool compute_sigma8;   // if true, then sigma8 will be computed, pars.add("output","mPk")

    bool use_SN;
    bool use_UNION;     // Union2.1 compilation
    bool use_SNLS;      // SNLS3 compilation
    bool use_SNLS_Mock; // SNLS3 mock sample
    bool use_JLA;       // Joint Lightcurve Analysis
    bool use_JLA_Mock;  // JLA mock sample
	bool use_WFIRST;	// WIFRST mock sample

//  ======
//  Priors
//  ======
    bool use_DDE_CPZ_prior;
    bool use_CMB_DistPrior; // CMB distance prior. no needed to do full CMB calculation.
    DDE_CPZ         *prior_dde_cpz;
    CMB_Dist        *prior_cmb_dist;

//	================
//	Data Structures:
//	================
    Data_Age		*data_age;
    Data_BAO        *data_bao;
    Data_HST        *data_HST;
    Data_Hubble     *data_Hz;
    Data_RSD        *data_rsd;

    SNe_UNION       *data_sne_union;
    SNe_SNLS        *data_sne_snls;
    SNe_SNLS_Mock   *data_sne_snls_mock;
    SNe_JLA         *data_sne_jla;
    SNe_JLA_Mock    *data_sne_jla_mock;
	SNe_WFIRST		*data_sne_wfirst;

#if defined(_USE_WMAP7_)
	Data_WMAP7		*data_wmap7;
#elif defined(_USE_WMAP9_)
    Data_WMAP9		*data_wmap9;
#endif
    Data_Planck2015 *data_plk2015;

    DataList();
    ~DataList();
    void Init( std::string& paramfile );
    void PrintParamList();

    void Get_MCMC_Params( std::vector<std::string>& mcmc_params );
    void Get_Derived_Params( std::vector<std::string>& derived_params );
};

#endif
