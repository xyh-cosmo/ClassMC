#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "ClassEngine.hh"
#include "CosmoTheory.hpp"
#include "CMB.hpp"
#include "class.h"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

Data_CMB_Dist::Data_CMB_Dist(){
    // do nothing
    use_Hu_fitting = true;
}

Data_CMB_Dist::~Data_CMB_Dist(){
    // do nothing
}


void Data_CMB_Dist::Init( string& CMB_dist_prior_dataset ){

    data_info.GetInfo(CMB_dist_prior_dataset);

    string dist_prior_file = "";
    string inv_covmat_file = "";
	string dist_prior_std_file = ""; // used only for Planck distance prior
	format = -1;	// 0--WMAP; 1--Planck

	if( Read::Has_Key_in_File(CMB_dist_prior_dataset,"format") ){
		format = Read::Read_Int_from_File(CMB_dist_prior_dataset,"format");
		if( format != 0 && format != 1 ){
			MPI_ClassMC_ERROR("format can only take 0 or 1 !");
		}
	}
	else{
		MPI_ClassMC_ERROR("the key \'format\' can not be found in file: "+CMB_dist_prior_dataset);
	}

    if( Read::Has_Key_in_File(CMB_dist_prior_dataset,"use_Hu_fitting") )
        use_Hu_fitting = Read::Read_Bool_from_File(CMB_dist_prior_dataset,"use_Hu_fitting");
    else{
        MPI_ClassMC_WARNING("use_Hu_fitting is not found, so use default value \'true\'.");
        // exit(0);
    }

    if( Read::Has_Key_in_File(CMB_dist_prior_dataset,"dist_prior_file") )
        dist_prior_file = Read::Read_String_from_File(CMB_dist_prior_dataset,"dist_prior_file");
    else
        MPI_ClassMC_ERROR("dist_prior_file not set in: "+CMB_dist_prior_dataset);


    if( Read::Has_Key_in_File(CMB_dist_prior_dataset,"inv_covmat_file") )
        inv_covmat_file = Read::Read_String_from_File(CMB_dist_prior_dataset,"inv_covmat_file");
    else
        MPI_ClassMC_ERROR("inv_covmat_file not set in: "+CMB_dist_prior_dataset);

	if( format == 0 ){
    	distance_prior  = arma::zeros(1,3);
    	covmat_inv      = arma::zeros(3,3);

    	distance_prior.load(dist_prior_file,arma::auto_detect); // loading CMB distance prior {lA,R,z_rec}
    	covmat_inv.load(inv_covmat_file,arma::auto_detect);     // loading covariance matrix
	}
	else if( format == 1 ){
//		cout << "loading Planck distance prior data ...\n";
		arma::rowvec tmp = arma::zeros(1,3);
		string std_file = Read::Read_String_from_File(CMB_dist_prior_dataset,"dist_prior_std_file");
		tmp.load(std_file,arma::auto_detect);
		std_R 	= tmp(0);
		std_lA	= tmp(1);
		std_Obh2= tmp(2);

		distance_prior_plk = arma::zeros(1,3);
		covmat_inv_plk     = arma::zeros(3,3);

    //  read Planck distance prior data
		distance_prior_plk.load(dist_prior_file,arma::auto_detect);

    //  read the normalized covariance matrix
		covmat_inv_plk.load(inv_covmat_file,arma::auto_detect);

    //  convert to the original covariance matrix
        covmat_inv_plk(0,0) *= std_R*std_R;
        covmat_inv_plk(0,1) *= std_R*std_lA;
        covmat_inv_plk(0,2) *= std_R*std_Obh2;

		covmat_inv_plk(1,0) *= std_lA*std_R;
		covmat_inv_plk(1,1) *= std_lA*std_lA;
		covmat_inv_plk(1,2) *= std_lA*std_Obh2;

		covmat_inv_plk(2,0) *= std_Obh2*std_R;
		covmat_inv_plk(2,1) *= std_Obh2*std_lA;
		covmat_inv_plk(2,2) *= std_Obh2*std_Obh2;

        covmat_inv_plk = covmat_inv_plk.i();
	}

}


double Prior_CMB_Dist(  imcmc_double&   param,
                        double&         lndet,
                        double&         chisq,
                        void*           model,
                        void*           data,
                        istate&         state ){

    chisq = lndet = 0.0;
    CosmoTheory *cosmo = static_cast<CosmoTheory*>(model);
    CMB_Dist    *cmbdist = static_cast<CMB_Dist*>(data);

    double lA, R, z_rec; // this format matches WMAP7/9 distance prior data

    double H0   = cosmo->engine->get_H0();
    double Om0  = cosmo->engine->get_Omega_m();
    double ra_rec, rs_rec;

    if( cmbdist->use_Hu_fitting ){
        z_rec  = cosmo->engine->z_rec_Hu();
        ra_rec = cosmo->engine->get_Da(z_rec)*(1.+z_rec)/cosmo->engine->get_a_today();
        rs_rec = cosmo->engine->get_rs(z_rec);
    }
    else{
        z_rec  = cosmo->engine->z_rec();
        ra_rec = cosmo->engine->ra_rec();
        rs_rec = cosmo->engine->rs_rec();
    }

    lA = _PI_*ra_rec/rs_rec;
    R  = sqrt(Om0)*H0*ra_rec/(_c_*1e-3);

    arma::rowvec residual = arma::zeros(1,3);

    switch (cmbdist->format){
        case 0:
        // for WMAP distance prior
            residual(0) = lA;
            residual(1) = R;
            residual(2)	= z_rec;
            residual = residual - cmbdist->distance_prior;
            chisq = arma::as_scalar( residual * cmbdist->covmat_inv * residual.t() );
            break;
        case 1:
        // for Planck distance prior
            residual(0) = R;
            residual(1) = lA;
            residual(2)	= cosmo->engine->get_Omega_b()*(H0/100)*(H0/100);
            residual = residual - cmbdist->distance_prior_plk;
            chisq = arma::as_scalar( residual * cmbdist->covmat_inv_plk * residual.t() );
            break;
        default:
            break;
    }

    return lndet - 0.5*chisq;
}
