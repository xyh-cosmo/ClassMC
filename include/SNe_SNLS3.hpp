//  ========================================================================
//  Author: Youhua Xu
//  Date:   Jan-7, 2017
//
//  Supernovae data structures and nuisance parameters
//
//  JLA is defined in another independent source file
//  ========================================================================

#ifndef __SNE_SNLS3__
#define __SNE_SNLS3__

#include <string>
#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>

#include "ParamList.hpp"


struct SNe_SNLS{

    DataInfo    data_info;

    int         sne_num;            //  472

    double 	    pecz;               // = 0.001 --> 300 km/sec

    bool        use_four_disp;      // if true, then use four dispersion for each subset
    double      intrinsicdisp;      //  0.13
    double      intrinsicdisp0;     //  0.0675
    double      intrinsicdisp1;     //  0.1133
    double      intrinsicdisp2;     //  0.0815
    double      intrinsicdisp3;     //  0.0989

    bool        twoscriptmfit;      //  default is set to false
    double      scriptmcut;         //  10
    double      scriptMB1, scriptMB2;   // used only if twoscriptmfit == True

    bool        marg_scriptm;   // default true, to be consistent with SNLS3 papers
    arma::vec   K1;
    arma::vec   K2;

    SNe_SNLS(); // do nothing
    ~SNe_SNLS(); // do nothing

//  lightcurvae data
    std::vector<std::string>  sne_name;
    std::vector<int>          set;

    arma::vec   zcmb, zhel, dz;
    arma::vec   mb, dmb;
    arma::vec   s, ds;
    arma::vec   c, dc;
    arma::vec   var3, dvar3;
    arma::vec   cov_m_s;
    arma::vec   cov_m_c;
    arma::vec   cov_s_c;

//  covariance matrix,
    arma::mat   cov_mB_mB;
    arma::mat   cov_mB_alpha;
    arma::mat   cov_mB_beta;
    arma::mat   cov_alpha_alpha;
    arma::mat   cov_alpha_beta;
    arma::mat   cov_beta_beta;

//  store pre-computed variable
    arma::vec   pre_vars;

//  every time alpha, beta, scriptMB been changed, the cov_tot will be updated
    double      alpha, beta, scriptMB;
    arma::mat   cov_tot, icov_tot;
    void        ReadCov( arma::mat& covmat, std::string& covmat_file );
    void        UpdateCov();
    void        UpdateNuisance( imcmc::imcmc_double par );
    void        SaveTotalCov( std::string tot_cov_fname );

    arma::vec   get_zcmb();
    arma::vec   get_zhel();
    double      max_z();

    void        Init( std::string& dataset );
};


struct SNe_SNLS_Mock{
    DataInfo    data_info;
    int         sne_num;
    bool        has_err;
    bool        use_full_cov;
	bool		invert_err; // used to check randomness of the direction of w transition
    arma::vec   z;
    arma::vec   mb, dmb;
    arma::mat   covmat_inv;
    SNe_SNLS_Mock();
    ~SNe_SNLS_Mock();
    void Init( std::string& dataset );
};


double Likelihood_SNe_SNLS( imcmc::imcmc_double&    param,
                            double&                 lndet,
                            double&                 chisq,
                            void*                   model,
                            void*                   data,
                            imcmc::istate&          state );

double Likelihood_SNe_SNLS_Mock(imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state );

#endif  //  __SNE_SNLS3__
