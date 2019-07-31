//	=================================================================
//	Wrapper to wmap likelihood functions written in damm fortran !!!
//	=================================================================

#ifndef __WMAP9__
#define __WMAP9__

#include <iostream>
#include <string>
#include <imcmc/imcmc.hpp>

#include "CosmoTheory.hpp"
#include "ParamList.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

extern "C"{
//	wmap options:
#ifdef USE_INTEL
void wmap_options_mp_wmap_print_options_(void);
void wmap_options_mp_wmap_set_options_(void);
#elif USE_GFORTRAN
void __wmap_options_MOD_wmap_print_options(void);
void __wmap_options_MOD_wmap_set_options(void);
#endif

//	wmap9 util:

#ifdef USE_INTEL
void wmap_util_mp_get_free_lun_(int *lun);
void wmap_util_mp_wmap_likelihood_error_init_(void);
void wmap_util_mp_wmap_likelihood_error_(int *msg, int *code);
void wmap_util_mp_wmap_likelihood_warning_(int *msg,int *code);
void wmap_util_mp_wmap_likelihood_error_report_(void);
#elif USE_GFORTRAN
void __wmap_util_MOD_get_free_lun(int *lun);
void __wmap_util_MOD_wmap_likelihood_error_init(void);
void __wmap_util_MOD_wmap_likelihood_error(int *msg, int *code);
void __wmap_util_MOD_wmap_likelihood_warning(int *msg, int *code);
void __wmap_util_MOD_wmap_likelihood_error_report(void);
#endif

//	wmap9 like:
#ifdef USE_INTEL
	void wmap_likelihood_9yr_mp_wmap_likelihood_init_(void);
	void wmap_likelihood_9yr_mp_wmap_likelihood_dof_(int *tt_npix, int *teeebb_npix);
	#if defined(USE_LOWELL_TBEB)
	void wmap_likelihood_9yr_mp_wmap_likelihood_compute_(double cltt[],double clte[],double cltb[],double clee[],double cleb[],double clbb[],double *like);
	#elif defined(USE_HIGHELL_TB)
	void wmap_likelihood_9yr_mp_wmap_likelihood_compute_(double cltt[],double clte[],double cltb[],double clee[],double cleb[],double clbb[],double *like);
	#else
	void wmap_likelihood_9yr_mp_wmap_likelihood_compute_(double cltt[],double clte[],double clee[],double clbb[],double *like);
	#endif
#elif USE_GFORTRAN
	void __wmap_likelihood_9yr_MOD_wmap_likelihood_init(void);
    void __wmap_likelihood_9yr_MOD_wmap_likelihood_dof(int *tt_npix, int *teeebb_npix);
    #if defined(USE_LOWELL_TBEB)
    void __wmap_likelihood_9yr_MOD_wmap_likelihood_compute(double cltt[],double clte[],double cltb[],double clee[],double cleb[],double clbb[],double *like);
    #elif defined(USE_HIGHELL_TB)
    void __wmap_likelihood_9yr_MOD_wmap_likelihood_compute(double cltt[],double clte[],double cltb[],double clee[],double cleb[],double clbb[],double *like);
	#else
	void __wmap_likelihood_9yr_MOD_wmap_likelihood_compute(double cltt[],double clte[],double clee[],double clbb[],double *like);
    #endif
#endif
}

struct Data_WMAP9{

	DataInfo data_info;

	Data_WMAP9();
	~Data_WMAP9();


	void Init( string& dataset );

    int lmax;   // 1200
	int tt_npix;
	int teeebb_npix;

};


double Likelihood_WMAP9(imcmc_double	&full_params,
						double			&lndet,
						double			&chisq,
						void			*model,
						void			*data,
						istate			&state );

#endif
