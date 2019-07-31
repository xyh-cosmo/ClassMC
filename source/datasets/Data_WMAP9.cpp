#include "WMAP9.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

Data_WMAP9::Data_WMAP9(){
    lmax = 1200;
}


Data_WMAP9::~Data_WMAP9(){
//	do nothing
}


void Data_WMAP9::Init( string& dataset ){
#ifdef USE_INTEL
	wmap_likelihood_9yr_mp_wmap_likelihood_init_();
	wmap_likelihood_9yr_mp_wmap_likelihood_dof_(&tt_npix,&teeebb_npix);
//    wmap_options_mp_wmap_print_options();
#elif USE_GFORTRAN
	__wmap_likelihood_9yr_MOD_wmap_likelihood_init();
	__wmap_likelihood_9yr_MOD_wmap_likelihood_dof(&tt_npix,&teeebb_npix);
 //   __wmap_options_MOD_wmap_print_options();
#endif

}
