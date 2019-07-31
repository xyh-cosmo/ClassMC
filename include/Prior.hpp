/*
	All priors play the same role of likelihoods do, so the prior functions
	defined here will incorporated into Likelihood(**)
*/

#ifndef __PRIOR__
#define __PRIOR__

#include <armadillo>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "class.h"
#include "Misc.hpp"
#include "ParamList.hpp"

//	====================================================================
//	Prior for reconstructing the Equation of State of Dark Energy : EoS
//	====================================================================

struct DDE_CPZ{

	//	Continuity prior with a const w^{fid} marginalized is called: Marginalized-Prior (MP)
	//	Continuity prior with w^{fid} approximated with local average is called: Smoothed-Prior (SP)

	DataInfo	data_info;

	bool 		floating;		//	if true, use local average to approximate w^{fid}

	int 		nw;             //  number of w bins or interpolation points
	int 		dde_spacing;	//	0 : evenly spaced in [0, zmax]
								//	1 : evenly spaced in [amin, 1]

	int 		smooth_bin_num; //  number of bins used to get (local-) smoothed w_fiducial

	std::string used_covmat_file, used_icovmat_file;	//	name of the output covmat and icovmat
	std::string	input_covmat_file;

	double 		*z, *a;
	double 		sigmam, zc;	    //	this is to be set from paramfile, sigmam_cpz and zc_cpz
    double 		sigmafid;       //  assumed fiducial variance of w^fid

	double 		sigma_mean_w;	//	\sigma_{\bar{w}}
	double 		ac;				//	correlation length in scale factor

	double 		dde_z_max;		//	max bin redshift
	double 		dde_a_min;		//	min redshifti

	double		zs, as;			//	smoothing scales used in Smoothed-Prior

	arma::mat   I;
	arma::mat	S;				//	used when floating=true.
	arma::mat 	C;	        	//	covariance matrix
	arma::mat 	iC;	        	//	inverse of the covariance matrix
	arma::mat	IS;				//	I-S

//	The following are added to perform w(a) reconstructions under LCDM assumption.
	// bool		use_lcdm_as_fiducial;
	arma::mat	C_original;		//	original covariance matrix proposed by CPZ
	arma::mat	iC_original;

	void Init( std::string& prior_settings );

	bool save_invcov;

    DDE_CPZ();
	~DDE_CPZ();

};

double Prior_DDE_CPZ(	imcmc::imcmc_double&    param,
                        double&                 lndet,
                        double&                 chisq,
                        void*                   model,
                        void*                   data,
                        imcmc::istate&          state );

#endif //__PRIOR__
