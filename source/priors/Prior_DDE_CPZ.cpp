//	Prior for reconstructing the Equation of State of Dark Energy
//	ref:
//
//  @Jan-19-2017: This file is mainly copied from another projetc of mine, xcos.
//  Some improvemnet were made to make the code more readable.
//
//	@May-17-2019: Now the covariance matrix is loaded from a pre-computed data file

#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "Prior.hpp"

using namespace imcmc;
using namespace imcmc::parser;
using namespace std;

DDE_CPZ::DDE_CPZ(){
    floating = false;
    z = NULL;
    a = NULL;
    save_invcov = false; // this is for debug, set to false by default
}

DDE_CPZ::~DDE_CPZ() {

    if( dde_spacing == 0 )
        delete[] a;
    else if( dde_spacing == 1 )
        delete[] z;

}

void DDE_CPZ::Init( std::string& prior_settings ) {

    data_info.GetInfo(prior_settings);

    floating		= Read::Read_Bool_from_File(prior_settings, "dde_wfid_floating");
    nw 				= Read::Read_Int_from_File(prior_settings, "eos_bin_num");

    dde_z_max 		= Read::Read_Double_from_File(prior_settings, "tabulated_w_z_max");
    dde_a_min		= 1./(1.+dde_z_max);

    sigmafid 		= Read::Read_Double_from_File(prior_settings, "dde_cpz_prior_sigmafid");

    save_invcov     = Read::Read_Bool_from_File(prior_settings, "save_invcov");

//	====================================
//  load pre-computed covariance matrix
    string cov_file = Read::Read_String_from_File(prior_settings, "pre-computed_covmat");

//	====================================
//  load pre-determined bin centers
    string abin_file= Read::Read_String_from_File(prior_settings, "abin_centers");

    arma::vec abin_centers;
    abin_centers.load(abin_file,arma::raw_ascii);

//	============================================================
//	@ Nov-8, 2016, update 'sigma_mean_w' option names ...
//	sigma_mean_w is the same for either evenly spaced in z or a
//	============================================================
    sigma_mean_w	= Read::Read_Double_from_File(prior_settings, "dde_cpz_prior_sigma_mean_w");

    C	= arma::zeros(nw, nw);
    iC	= arma::zeros(nw, nw);

    C.load(cov_file, arma::raw_ascii);

    //  We choose sigma_mean_w=0.04 as a reference , so if one wants to change sigma_mean_w,
    //  he just needs to rescale the loaded covariance matrix without re-compute
    //  the covariance matrix
    C = C * pow(sigma_mean_w/0.04,2);

    iC	= C.i();	//	evaluate the inverse of the covariance matrix

//  make a copy of the original covariance matrix without local-average
    C_original  = C;
    iC_original = iC;

    dde_spacing = Read::Read_Int_from_File(prior_settings,"dde_spacing");
    if( dde_spacing == 0 )
        a = new double[nw];
    else if( dde_spacing == 1 )
        z = new double[nw];

//	====================================================================
//	NOTE: zbin[] and abin[] defined below is DIFFERENT from that used in
//	eos approximation !!!
//	====================================================================
    double zbin[nw];
    double abin[nw];

    double Delta;
    double xi0;
    double xij, x_plus, x_minus, xbar;

    if( dde_spacing == 0 ) {		//	evenly spaced in a

        ac	= Read::Read_Double_from_File(prior_settings, "dde_cpz_prior_ac");

        //  try to read smoothing scale $a_s$ from "prior_settings"
		if( Read::Has_Key_in_File(prior_settings,"dde_cpz_prior_as") ){
			as	= Read::Read_Double_from_File(prior_settings, "dde_cpz_prior_as");
		}
		else{
			as	= ac;
		}

    //  Correct bin center of the last abin to make sure that it is correlated with the one next to it.
    // 	Added by XYH @20190516
        abin_centers[nw-1] = abin_centers[nw-2] - (as-1e-10);

        double da = (1.0-dde_a_min)/nw;
        for(int i=0; i<nw; ++i) {
            abin[i] = 1.0 - (i+0.5)*da; // center of each a-bin

            if( floating ){
                a[i] = abin_centers[i];
            }
        }

        Delta 	= da;
        xi0		= sigma_mean_w*sigma_mean_w*(1-dde_a_min)/_PI_/ac;
    }
    else if( dde_spacing == 1 ) {	//	evenly spaced in z

        zc 	= Read::Read_Double_from_File(prior_settings, "dde_cpz_prior_zc");

		if( Read::Has_Key_in_File(prior_settings,"dde_cpz_prior_zs") ) {
			zs	= Read::Read_Double_from_File(prior_settings, "dde_cpz_prior_zs");
		}
		else{
			zs	= zc;
		}

        double dz = dde_z_max/nw;
        for(int i=0; i<nw; ++i) {
            zbin[i] = (i+0.5)*dz;	// center of each a-bin

            if( floating )
                z[i] = zbin[i];
        }

        Delta 	= dz;
        xi0		= sigma_mean_w*sigma_mean_w*dde_z_max/_PI_/zc;
    }
    else {
        string err = "\n*** DDE_CPZ::Init() ==> unsupported z-spacing !!!";
        throw runtime_error(err);
    }

    S	= arma::zeros(nw, nw);  // smoothing matrix
    I	= arma::eye(nw,nw);

    if( floating ) {

        for(int i=0; i<nw; ++i) {
            int count=0;

            if( dde_spacing == 0 ) {
                for( int j=0; j<nw; ++j ) {
                    if( fabs(a[j] - a[i]) <= as + 1E-10 ) {
                        S(i,j) = 1.0;
                        ++count;
                    }
                }
            }
            else if( dde_spacing == 1 ) {
                for( int j=0; j<nw; ++j ) {
                    if( fabs(z[j] - z[i]) <= zs + 1E-10 ) {
                        S(i,j) = 1.0;
                        ++count;
                    }
                }
            }

            for( int j=0; j<nw; ++j )
                S(i,j) = S(i,j)/count;
        }



        IS = I-S;
        iC = IS.t() * iC * IS;
    }

//    S.save("EoS_smooth_matrix.txt", arma::raw_ascii);
//    iC.save("EoS_inv_covmat_20190701.txt",arma::raw_ascii);

    if( save_invcov ){
        if( MPI::COMM_WORLD.Get_rank() == 0 ){
            cout << "Saving inverse of the CPZ covariance matrix \n";
            iC.save("DDE_CPZ_prior_invcov_uncorrected.txt",arma::raw_ascii);
        }
    }

//    exit(0);
}

double Prior_DDE_CPZ(  	imcmc_double&   param,
                        double&         lndet,
                        double&         chisq,
                        void*           model,
                        void*           data,
                        istate&         state ) {

    state.this_like_is_ok=true;

    lndet = chisq = 0;

    DDE_CPZ *cpz = static_cast<DDE_CPZ*>(data);

    arma::rowvec row_w = arma::zeros<arma::rowvec>(cpz->nw);
    arma::colvec col_w = arma::zeros<arma::colvec>(cpz->nw);

//  those extra high-z w_i will not be constrained by CPZ prior.

    for(int i=0; i<cpz->nw; ++i) {
        string wname= "DDE_w"+Read::IntToString(i);
        row_w[i] = param[wname];
        col_w[i] = row_w[i];
    }

    chisq = arma::as_scalar( row_w * cpz->iC * col_w );

    return -lndet - 0.5*chisq;
}
