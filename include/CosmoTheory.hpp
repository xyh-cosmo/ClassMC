#ifndef __COSMOTHEORY__
#define __COSMOTHEORY__

#include <imcmc/parser++.hpp>
#include <imcmc/ensemble.hpp>

#include "Engine.hh"
#include "ClassEngine.hh"
#include "Misc.hpp"

using namespace imcmc;
using namespace imcmc::parser;

#if defined(_USE_WMAP7_)
    struct Data_WMAP7;
#elif defined(_USE_WMAP9_)
    struct Data_WMAP9;
#endif

struct Data_Planck2015;
struct DataList;

struct CosmoTheory{

//  ============================================================================
//  pars is used to tell CLASS what parameters will be used, but only in the
//  initialization of engine. After the engine's initialization, filecontent FC
//  will have recorded all the necessary names of the parameters needed by CLASS.
//  Then during MCMC sampling, we uss imcmc_double& full_param to update the
//  filecontent FC.
    ClassParams     pars;
    ClassEngine     *engine;

    bool TT, EE, BB, TE, TB, EB;
    bool do_CMB;    // if no CMB calculation or CMB data is not used, set this to false

//  this flag is used to control the calculation of thermoal history. E.g., when using
//  Hu-fitting formula and CMB distance prior, there is no need to compute the whole
//  recombination history
    bool do_thermo;

//  Cl at high l's are suppressed by weak lensing effect.
    bool lensed_CMB;

//  what do you want as the outputs?
    std::string  output;

//  ============================================================================
//  mcmc_params makes a record of all the MCMC parameters, and will tell imcmc
//  what to read from a *ini file.
    std::vector<std::string> mcmc_params;

//  derived parameters. Always be careful with derived parameters !!!!!
    bool has_derived_params;
    std::vector<std::string> derived_params;
//  read from 'derived_cosmo_params'
    void add_derived_params( string& paramfile );
    void update_derived_params( imcmc_double& full_param );

//  ============================================================================
//  Member functions
    void Init( std::string& paramfile );
    void Init( Data_Planck2015* clik );
#if defined(_USE_WMAP7_)
	void Init( Data_WMAP7 *wmap7 );
#elif defined(_USE_WMAP9_)
    void Init( Data_WMAP9 *wmap9 );
#endif
    void Init( DataList& datalist );
    bool update( imcmc_double& full_param );

//  ============================================================================
//  constructor and destructor
    CosmoTheory();
    ~CosmoTheory();

//  ============================================================================
//  Some options controlling CLASS

//  opt_As controls the parameterization of As:
    int opt_As; //  0 = A_s; 1 = 10^9A_s; 2 = ln(10^10A_s)

//  use_sBBN controls wether use BBN model
//  if true, YHe will be treated as derived parameter (default: true)
    bool use_sBBN;

//  if true, low-z EoS(z) of dark energy is approximated by interpolation
//  over {zi,wi}
    bool use_DDE;	// this option is needed by CLASS

//  ============================================================================
//  MCMC parameters: parameters that denote the same quanties which are used both
//  in class and MCMC sampling may be different, so we need to pair vars in class
//  with those used to comunicate parameters between imcmc and likelihoods.
    std::map<std::string,std::string> class_mcmc_params;
    void Pair_mcmc_params( std::string class_mcmc_pname, std::string imcmc_pname, double value );
    void Print_pair_mcmc_params();
    void Get_Class_MCMC_Params( std::vector<std::string>& class_mcmc_prams );
    void Get_Derived_Params( std::vector<std::string>& derived_params );

};


struct CosmoTheoryTest{

    bool test_mu;   // distance modulus
    bool test_BAO;  // BAO. 5 observables.
    bool test_Hz;   // H(z)
    bool test_weff; // effective EoS

    int output_size;
    double zmin,zmax;

    std::string test_root;

    void Init(std::string& paramfile);
    void RunTest(std::string& paramfile, CosmoTheory *cosmo);
};


#endif  //__COSMOTHEORY__
