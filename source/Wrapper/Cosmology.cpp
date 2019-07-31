//  =================================================================================
//  This source file contains the CMB calculator (based on CLASS) for plik
//
//  Note:
//  Knowing how model parameters are passed between Class and IMCMC is important!!!
//  1) 'pars' is of the type ClassParams, it will be used to create and initialize
//     file_content pointers, which will used by imcmc_double full_params to update
//     newly proposed MCMC sampling parameters.
//  2) Also note that not all parameters stored in pars (or file_content pointer) are
//     MCMC parameters, however, we still make a copy of them and stroe them into full_params.
//     The reason is that updating the file_content pointer can be easily done by
//     re-copying all parameters stored in full_params into the file_content pointer.
//  =================================================================================

#include "ClassMC.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

CosmoTheory::CosmoTheory() {
    engine      = NULL;
    lensed_CMB  = false;
    output      = "";
    use_sBBN    = false;
    use_DDE     = false;
    has_derived_params = false;
    do_CMB      = false;
    do_thermo   = false;
}

CosmoTheory::~CosmoTheory() {
    if( engine != NULL ) {
        delete engine;
        engine = NULL;
    }
}

void CosmoTheory::Pair_mcmc_params( string class_mcmc_pname,
                                    string imcmc_pname,
                                    double value ){
    if( class_mcmc_params.count(class_mcmc_pname) == 0 ){
        pars.add(class_mcmc_pname, value);
        mcmc_params.push_back(imcmc_pname);
        class_mcmc_params[class_mcmc_pname] = imcmc_pname;
    }
    else{
        string errmsg = class_mcmc_pname + " already existed in Class_MCMC_Params!";
		MPI_ClassMC_ERROR(errmsg);
    }
}

void CosmoTheory::Print_pair_mcmc_params(){
    map<string,string>::iterator it;
    for( it=class_mcmc_params.begin(); it!=class_mcmc_params.end(); ++it ){
        cout << "Class_MCMC_Params: " << it->first << " <--> " << it->second << endl;
    }
}

//  =================================================================================
//  if you choose different cosmological parametrizations, MODIFY this function !!!!!
//  =================================================================================

void CosmoTheory::Init( string& paramfile ) {

    class_mcmc_params.clear();
    mcmc_params.clear();

    if( Read::Has_Key_in_File(paramfile,"do_thermo") ){
        do_thermo = Read::Read_Bool_from_File(paramfile,"do_thermo");
    }

    double value;    //    we are only dealing with vars of "double"-type

    if( Read::Has_Value(paramfile, "100*theta_s", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "100*theta_s");
        Pair_mcmc_params("100*theta_s","100*theta_s",value);
    }
    else if( Read::Has_Value(paramfile, "H0", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "H0");
        Pair_mcmc_params("H0","H0",value);
    }
    else {
        MPI_ClassMC_ERROR("## you must provide either 100*theta_s or H0 !!!");
    }

    //  Baryons
    //  Note the conventions used here: 
    //  'omegab'   and 'Omega_b' ---> $\Omega_b$
    //  'omegabh2' and 'omega_b' ---> $\Omega_b h^2$

    if( Read::Has_Value(paramfile, "omegab", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "omegab");
        Pair_mcmc_params("Omega_b","omegab",value);
    }
    else if( Read::Has_Value(paramfile, "omegabh2", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "omegabh2");
        Pair_mcmc_params("omega_b","omegabh2",value);
    }
    else {
        MPI_ClassMC_ERROR("## you must provide either omegab or omegabh2 !!!");
    }

    //  Ultra-relativistic matter. (we consider only neutrinos)

    double N_ur = 3.046;
    if( Read::Has_Value(paramfile, "N_ur", "double") ) {
        N_ur = Read::Read_Double_from_File(paramfile, "N_ur");
        Pair_mcmc_params("N_ur","N_ur",N_ur);
    }
    else {
        pars.add("N_ur", N_ur);
    }

    //  Cold dark matter

    if( Read::Has_Value(paramfile, "omegac", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "omegac");
        Pair_mcmc_params("Omega_cdm","omegac",value);
    }
    else if( Read::Has_Value(paramfile, "omegach2", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "omegach2");
        Pair_mcmc_params("omega_cdm","omegach2",value);
    }
    else {
        MPI_ClassMC_ERROR("## you must provide either omegab or omegach2 !!!");
    }

    //  Curvature

    if( Read::Has_Value(paramfile, "omegak", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "omegak");
        Pair_mcmc_params("Omega_k","omegak",value);
    }

//    ==========================================================================
//    Massive neutrions should be carefully taken care of ...
    int N_ncdm = 0;
    // double m_ncdm;

    if( Read::Has_Value(paramfile, "N_ncdm", "int") ) {
        N_ncdm = Read::Read_Int_from_File(paramfile, "N_ncdm");
        Pair_mcmc_params("N_ncdm","N_ncdm",N_ncdm);
    }
    else {   //    set to zero, this is the default case
        pars.add("N_ncdm", 0);
    }

//  ==================================
//  This part needs improvements !!!!!
//  ==================================
    if( N_ncdm > 0 ) {   //    could be 1, 2 or 3

        // cout << "N_ncdm = " << N_ncdm << endl;

        if( N_ncdm == 1 ) {
            pars.add("m_ncdm", "0.06");
            // mcmc_params.push_back("m_ncdm_tot");    //    m_ncdm = m_ncdm_tot/N_ncdm
            cout << "# Your universe has 1 massive neutrino!\n";
        }
        else if( N_ncdm == 2 ) {
            pars.add("m_ncdm", "0.03, 0.03");
            // mcmc_params.push_back("m_ncdm_tot");    //    m_ncdm = m_ncdm_tot/N_ncdm
            cout << "# Your universe has 2 massive neutrinos!\n";
        }
        else if( N_ncdm == 3 ) {
            pars.add("m_ncdm", "0.02, 0.02, 0.02");
            // mcmc_params.push_back("m_ncdm_tot");    //    m_ncdm = m_ncdm_tot/N_ncdm
            cout << "# Your universe has 3 massive neutrinos!\n";
        }
        else {
            MPI_ClassMC_ERROR("***** invalid setting of m_ncdm");
        }
    }

    //    check total Neff [[ improvements are needed here ! ]]
    if( (N_ur + N_ncdm) > 4.5) {
        //  to be more specific, N_ncdm is only approximately equal to Neff_{ncdm}
        MPI_ClassMC_ERROR("## wired Neff(N_ur+N_ncdm) = " + str(N_ur + N_ncdm) );
    }

    //  ===================
    //  Dark Energy part:
    //  -------------------
    //  == case 0):
    //  LCDM corresponds to setting Omega_fld=0, since Omegal is to be infered
    //  from other parametes, so there is no need to explicitly specify Omegal.
    //  == case 1):
    //  if want to use w(a)=w0+wa(1-a), then unset Omega_fld, so that it will be
    //  inferred from other parameters; at the same time, Omega_lambda must be
    //  set to zero, so that class can corretly know how to model the dark energy.
    //  == case 2):
    //  use scalar field model for DE. set Omega_lambda and Omega_fld to zero.
    //  this part has not tested ...
    //  == case 3):
    //  use DDE approxiamtion. set use_DDE to 1 is ok.

    int DE_option = Read::Read_Int_from_File(paramfile,"DE_option");

    if( DE_option == 0 ) {
        pars.add("Omega_fld", 0.0);
    }
    else if( DE_option == 1 ) {

        pars.add("Omega_Lambda", 0.0);

        value = -1.0;
        if( Read::Has_Value(paramfile,"w0_fld", "double") ){
            value = Read::Read_Double_from_File(paramfile, "w0_fld");
            Pair_mcmc_params("w0_fld","w0_fld",value);
        }

        value = 0.0;
        if( Read::Has_Value(paramfile,"wa_fld", "double") ){
            value = Read::Read_Double_from_File(paramfile, "wa_fld");
            Pair_mcmc_params("wa_fld","wa_fld",value);
        }

        value = 1.0;
        if( Read::Has_Value(paramfile,"cs2_fld", "double") ) {
            value = Read::Read_Double_from_File(paramfile, "cs2_fld");
            Pair_mcmc_params("cs2_fld","cs2_fld",value);
        }
    }
    else if( DE_option == 2 ) {
        MPI_ClassMC_ERROR("==>> This part has not been implemented and test yet ... <<==");
    }
    else if( DE_option == 3 ) {

        int DDE_table_size = -1;
        DDE_table_size       = Read::Num_of_Value_for_Key(paramfile,"DDE_z");

        if( DDE_table_size <= 0) {
            MPI_ClassMC_ERROR("==> Error in getting size of DDE_table (DDE_z or DDE_w)");
        }

    //  =================
    //  add {w_i} for DDE
        for( int i=0; i<DDE_table_size; ++i ) {
            std::string DDE_wi = "DDE_w"+Read::IntToString(i);
            double DDE_w = Read::Read_Double_from_File(paramfile,DDE_wi);
            Pair_mcmc_params(DDE_wi,DDE_wi,DDE_w);
        }

    //  =======================
    //  add extra {w_i} for DDE
        int DDE_table_size_extra;
        if( Read::Has_Key_in_File(paramfile,"DDE_z_extra") )
            DDE_table_size_extra = Read::Num_of_Value_for_Key(paramfile,"DDE_z_extra");
        else
            DDE_table_size_extra = 0;

        for( int i=0; i<DDE_table_size_extra; ++i ){
            std::string DDE_w_extra_i = "DDE_w_extra_"+Read::IntToString(i);
            double DDE_w = Read::Read_Double_from_File(paramfile,DDE_w_extra_i);
            Pair_mcmc_params(DDE_w_extra_i,DDE_w_extra_i,DDE_w);
        }

    //  =================
    //  add {z_i} for DDE
        double *DDE_z_values = new double[DDE_table_size];
        Read::Read_Array_of_Double_from_File(paramfile,"DDE_z",DDE_z_values,DDE_table_size);

        std::string DDE_z="";
        for( int i=0; i<DDE_table_size; ++i ) {
            std::string temp = Read::DoubleToString(DDE_z_values[i],6);
            if( i != DDE_table_size-1 )
                DDE_z += temp+",";
            else
                DDE_z += temp;
        }

        delete[] DDE_z_values;

        pars.add("DDE_z",DDE_z);
        pars.add("use_DDE",1);

//#if defined(_DEBUG_WEFF_)
//        MPI_cout("running weff debug initialization part ...\n");
//        double *DDE_w_values = new double[DDE_table_size+DDE_table_size_extra];
//        Read::Read_Array_of_Double_from_File(paramfile,"DDE_w",DDE_w_values,DDE_table_size);

//        std::string DDE_w="";
//        for( int i=0; i<DDE_table_size; ++i ) {
//            std::string temp = Read::DoubleToString(DDE_w_values[i],6);
//            if( i != DDE_table_size-1 )
//                DDE_w += temp+",";
//            else
//                DDE_w += temp;
//        }

//        delete[] DDE_w_values;
//        pars.add("DDE_w",DDE_w);
//        pars.add("DDE_w_format",0); // this is for debug and DDE_w[] is read from DDE_w = ..., instead of from DDE_wi = ..
//#else
//        pars.add("DDE_w_format",1); // tell class to read DDE_wi instead of DDE_w, very important !!
//#endif

		pars.add("DDE_w_format",1); // tell class to read DDE_wi instead of DDE_w, very important !!

    //  =======================
    //  add extra {z_i} for DDE
        if( DDE_table_size_extra > 0 ){
            std::string DDE_z_extra="";

            double *DDE_z_extra_values = new double[DDE_table_size_extra];
            Read::Read_Array_of_Double_from_File(paramfile,"DDE_z_extra",
                                                    DDE_z_extra_values,
                                                    DDE_table_size_extra);

            for( int i=0; i<DDE_table_size_extra; ++i ) {
                std::string temp = Read::DoubleToString(DDE_z_extra_values[i],6);
                if( i != DDE_table_size_extra-1 )
                    DDE_z_extra += temp+",";
                else
                    DDE_z_extra += temp;
            }

            delete[] DDE_z_extra_values;
            pars.add("DDE_z_extra",DDE_z_extra);
        }

    //  =====================================================
    //  trying to find some flag that controlling GSL spline
        if( Read::Has_Key_in_File(paramfile,"spline_DDE_weff_size") ){
            int spline_DDE_weff_size;
            spline_DDE_weff_size = Read::Read_Int_from_File(paramfile,"spline_DDE_weff_size");
            pars.add("spline_DDE_weff_size",spline_DDE_weff_size);
        }

        int gsl_interp_method;
        if( Read::Has_Key_in_File(paramfile,"w_gsl_interp_method") ){
            gsl_interp_method = Read::Read_Int_from_File(paramfile,"w_gsl_interp_method");
            pars.add("w_gsl_interp_method",gsl_interp_method);
        }

        if( Read::Has_Key_in_File(paramfile,"weff_gsl_interp_method") ){
            gsl_interp_method = Read::Read_Int_from_File(paramfile,"weff_gsl_interp_method");
            pars.add("weff_gsl_interp_method",gsl_interp_method);
        }

        if( Read::Has_Key_in_File(paramfile,"romberg_weff_eps") ){
            double romberg_weff_eps;
            romberg_weff_eps = Read::Read_Double_from_File(paramfile,"romberg_weff_eps");
            pars.add("romberg_weff_eps",romberg_weff_eps);
        }

    //  configure the approximation method: interpolation (==0) or piecewise constant (==1)
        if( Read::Has_Key_in_File(paramfile,"DDE_approx_opt") ){
            int DDE_approx_opt = Read::Read_Int_from_File(paramfile,"DDE_approx_opt");
            pars.add("DDE_approx_opt",DDE_approx_opt);
        }
        else{
            pars.add("DDE_approx_opt",0); // we use interpolation method by default.
        }
        
        double DDE_z_max = Read::Read_Double_from_File(paramfile,"DDE_z_max");
        pars.add("DDE_z_max",DDE_z_max);
    }

//  =================================================
//  configure which As parameterization is to be used
//  =================================================  
    opt_As = -1;

    if( Read::Has_Value(paramfile, "A_s", "double") ) { //  opt_As=0
        value = Read::Read_Double_from_File(paramfile, "A_s");
        Pair_mcmc_params("A_s","A_s",value);
        opt_As = 0;
    }
    else if( Read::Has_Value(paramfile, "10^9A_s", "double") ) {   //  opt_As=1
        value = Read::Read_Double_from_File(paramfile, "10^9A_s");
        Pair_mcmc_params("A_s","10^9A_s",value*1.E-9);
        opt_As = 1;
    }
    else if( Read::Has_Value(paramfile, "ln(10^10A_s)", "double") ) {       //  opt_As=2
        value = Read::Read_Double_from_File(paramfile, "ln(10^10A_s)");
        Pair_mcmc_params("A_s","ln(10^10A_s)",exp(value)*1.E-10);
        opt_As = 2;
    }
    else {
        MPI_ClassMC_ERROR("### you must provide at least one of {A_s, 10^9A_s, ln(10^10A_s}");
    }

//  spectrual index of the primodial power spectrum
    if( Read::Has_Value(paramfile, "n_s", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "n_s");
        Pair_mcmc_params("n_s","n_s",value);
    }

//  running of n_s
    if( Read::Has_Value(paramfile, "alpha_s", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "alpha_s");
        Pair_mcmc_params("alpha_s","alpha_s",value);
    }

    if( Read::Has_Value(paramfile, "tau_reio", "double") ) {    // optical depth to reionization
        value = Read::Read_Double_from_File(paramfile, "tau_reio");
        Pair_mcmc_params("tau_reio","tau_reio",value);
    }
    else if( Read::Has_Value(paramfile, "z_reio", "double") ) { // redshift of reionization
        value = Read::Read_Double_from_File(paramfile, "z_reio");
        Pair_mcmc_params("z_reio","z_reio",value);
    }

    if( Read::Has_Key_in_File(paramfile, "use_sBBN") ) {
        use_sBBN = Read::Read_Bool_from_File(paramfile,"use_sBBN");
    }

    if( use_sBBN ) { // YHe derived from sBBN is denoted as YHe_bbn
        pars.add("YHe","BBN");
        cout << "==> assuming standard BBN.\n";
        // exit(0);
    }
    else {
        if( Read::Has_Value(paramfile, "YHe", "double") ) {
            value = Read::Read_Double_from_File(paramfile, "YHe");
            Pair_mcmc_params("YHe","YHe",value);
            cout << "==> YHe is treated as a free parameter.\n";
        }
        else {
            pars.add("YHe",0.255);
            cout << "==> YHe is set to default value: 0.255\n";
        }
    }

//  ============================================================================
//  the following will not be sampled, but you can also add/remove some of them
//  NOTE: All NON-MCMC parameters MUST be placed at the end of pars, this is
//  extremly important !!!
//  ============================================================================


//  reionization_exponent = 1.5
    if( Read::Has_Value(paramfile, "reionization_exponent", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "reionization_exponent");
        pars.add("reionization_exponent", value);    //    no mcmc parameter
    }
    else {
        pars.add("reionization_exponent", 1.5);
    }

//  reionization_width = 1.5
    if( Read::Has_Value(paramfile, "reionization_width", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "reionization_width");
        pars.add("reionization_width", value);    //    no mcmc parameter
    }
    else {
        pars.add("reionization_width", 1.5);
    }

//  helium_fullreio_redshift = 3.5
    if( Read::Has_Value(paramfile, "helium_fullreio_redshift", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "helium_fullreio_redshift");
        pars.add("helium_fullreio_redshift", value);    //    no mcmc parameter
    }
    else {
        pars.add("helium_fullreio_redshift", 3.5);
    }

//  helium_fullreio_width = 0.5
    if( Read::Has_Value(paramfile, "helium_fullreio_width", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "helium_fullreio_width");
        pars.add("helium_fullreio_width", value);    //    no mcmc parameter
    }
    else {
        pars.add("helium_fullreio_width", 0.5);
    }

    if( Read::Has_Value(paramfile, "k_pivot", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "k_pivot");
        pars.add("k_pivot",value);
    }
    else {
        pars.add("k_pivot",0.05);
    }

    // T_cmb
    if( Read::Has_Value(paramfile, "T_cmb", "double") ) {
        value = Read::Read_Double_from_File(paramfile, "T_cmb");
        pars.add("T_cmb",value);
    }
    else {
        pars.add("T_cmb",2.7255);
    }

    if( Read::Has_Value(paramfile, "class_root", "string") ) {
        string class_root = Read::Read_String_from_File(paramfile, "class_root");
        pars.add("root", class_root);
    }

    if( Read::Read_Bool_from_File(paramfile, "write_class_params") == true ) {
        pars.add("write parameters", "yes");
    }

    //  ==============================
    //  add derived parameters if any
    //  ==============================
    has_derived_params = false;
    add_derived_params( paramfile ); // add derived parameter from key: derived_cosmo_params

    if( use_sBBN ) {
        bool has_YHe_bbn = false;

        for( size_t i=0; i<derived_params.size(); ++i ) {
            if( derived_params[i] == "YHe_bbn" ) {
                has_YHe_bbn = true;
                break;
            }
        }

        if( !has_YHe_bbn)
            derived_params.push_back("YHe_bbn");
    }

    // exit(0);

    //  =====================
    //  recombination option
    //  =====================
    if( Read::Has_Key_in_File(paramfile,"recombination") ) {
        string rec = Read::Read_String_from_File(paramfile,"recombination");
        if( rec != "" ) {
            pars.add("recombination", rec);
            cout << "CosmoTheory::init() ==> recombination option: " << rec << endl;
        }
        else {
            cout << "CosmoTheory::init() ==> key recombination has no value, so use default option\n"
                 << "RECFAST to compute recombination (see class document for more details).\n";
        }
    }
    else {
        cout << "CosmoTheory::init() ==> recombination is not found, so we will use default option\n"
             << "RECFAST to compute recombination (see class document for more details).\n";
    }
//    exit(0);

//  added @ Dec-18-2016:
//  NOTE: if lensing effects of Cls are not accounted for, the constrained cosmological
//  parameters are biased when using PLK.
    lensed_CMB = Read::Read_Bool_from_File(paramfile,"lensed_CMB");
}

//  determine what to be computed...
void CosmoTheory::Init( Data_Planck2015 *clik ) {

    bool has_tCl = false;
    bool has_pCl = false;
    bool has_lCl = false;

    output = "";        // default compute temperature & polarization
    int l_max_scalars = 2; // minmum value

    if( clik != NULL ) { // if plik is not used, Planck data pointer inside DataList will not be allocated

        clik->get_Cl_tasks(has_tCl,has_pCl,has_lCl);

        if( clik->has_highl ) {
            for( int i=0; i<6; ++i ) {
                if( clik->lmax_highl[i] > l_max_scalars )
                    l_max_scalars = clik->lmax_highl[i];
            }
        }

        if( clik->has_lowP ) {
            for( int i=0; i<6; ++i ) {
                if( clik->lmax_lowP[i] > l_max_scalars )
                    l_max_scalars = clik->lmax_lowP[i];
            }
        }

        if( clik->has_lens ) {

            for( int i=0; i<7; ++i ) {
                if( clik->lmax_lens[i] > l_max_scalars )
                    l_max_scalars = clik->lmax_lens[i];
            }
        }


        if( has_tCl )
            output += "tCl,";

        if( has_pCl )
            output += "pCl,";

        if( has_lCl || lensed_CMB ) {
            if( has_tCl || has_pCl )
                output += "lCl";
            pars.add("lensing", true); //note boolean
        }

		output += ",mPk";	// mPk is needed to compute sigma8

        pars.add("output",output.c_str());
        pars.add("l_max_scalars",l_max_scalars);
    }

}

#if defined(_USE_WMAP7_)
void CosmoTheory::Init( Data_WMAP7 *wmap7 ){
//  in CosmoMC, only 4 Cls are used: TT, TE, EE, BB
    output = "tCl,pCl,lCl,mPk";
    pars.add("output",output.c_str());
    pars.add("l_max_scalars",wmap7->lmax);
    pars.add("lensing", true); // here we always compute lensed power spectrum
}
#elif defined(_USE_WMAP9_)
void CosmoTheory::Init( Data_WMAP9 *wmap9 ){
//  in CosmoMC, only 4 Cls are used: TT, TE, EE, BB
    output = "tCl,pCl,lCl,mPk";
    pars.add("output",output.c_str());
    pars.add("l_max_scalars",wmap9->lmax);
    pars.add("lensing", true); // here we always compute lensed power spectrum
}
#endif

void CosmoTheory::Init( DataList& datalist ) {

	if( datalist.use_CMB ){

		if( datalist.use_PLK )
    		Init(datalist.data_plk2015);
#if defined(_USE_WMAP7_)
		else if( datalist.use_WMAP7 )
			Init(datalist.data_wmap7);
#elif defined(_USE_WMAP9_)
		else if( datalist.use_WMAP9 )
			Init(datalist.data_wmap9);
#endif
		pars.add("format","CAMB");
		do_CMB = true;
        do_thermo = true;
	}
	else{
        if( datalist.use_RSD || datalist.compute_sigma8 ){
            cout << "using redshift distortion measurements\n";
//            pars.add("output","tCl,pCl,lCl,mPk"); // calculation of sigma8 needs mPk  tCl,pCl,lCl,mPk
            pars.add("output","mPk"); // calculation of sigma8 needs mPk  tCl,pCl,lCl,mPk
            pars.add("l_max_scalars",2); // changing this number has very small effect to the sigma8 results.
//            pars.add("l_max_scalars",1200);
//            pars.add("l_max_tensors",500);
//            pars.add("l_max_lss",600);
//            pars.add("z_pk","0.01,0.02,0.03,0.04,0.05, 0.067, 0.08, 0.10");
            pars.add("z_max_pk",3.0); // this is fixed temporarily
            pars.add("P_k_max_h/Mpc", 1.);
//            do_CMB = false;
            do_CMB = true;  // this flag must be turned to assure correct sigma8 calculations.
            do_thermo = true;
            
//            printf("@@@@@  using compute_sigma8 @@@@@\n");
        }
		else{
        // in this case, sigma8 can not be computed (added @Aug-6-2017)
            pars.add("output","tCl");
            pars.add("l_max_scalars",2);
            do_CMB = false;
            do_thermo = false;
        }
	}

//  theromal thistory is needed for computing CMB & determing z_rec, z_drag.
    if( do_CMB == false ){

    //  if we use class to get rs_zrec, then need to set do_thermo to true.
        if( datalist.use_CMB_DistPrior && (!datalist.prior_cmb_dist->use_Hu_fitting) ){
            do_thermo = true;
        }

    //  if we use class to get rs_zdrag, then need to set do_thermo to true
        if( datalist.use_BAO && (!datalist.data_bao->use_Hu_fitting) ){
            do_thermo = true;
        }
    }

//  ====================
//  allocate ClassEngine
//  ====================
    engine = new ClassEngine(pars);
    engine->do_CMB      = do_CMB;
    engine->do_thermo   = do_thermo;  //  if use BAO, then we need to know the thermo history
    engine->use_sBBN    = this->use_sBBN;

//  ==================================================================
//  Now add cosmological parameters into datalist.ParList->MCMC_Params.
//  An important reason I put this step here is that, if ClassEgine is
//  successfully initialized, then the set of cosmological parameters
//  MUST form a complete set.

    MPI_cout("CosmoTheory::Init() ==> adding cosmological parameters into MCMC params ...");

    if( mcmc_params.size() <= 0 ){
        MPI_ClassMC_ERROR("It's wired that your mcmc_params is empty !");
    }

    for( size_t i=0; i<mcmc_params.size(); ++i ){
        datalist.ParList->Add_MCMC_Param(mcmc_params[i]);
        // cout << "adding " << mcmc_params[i] << endl;
    }

    MPI_cout("CosmoTheory::Init() ==> adding derived cosmological parameters into derived params ...");
    for( size_t i=0; i<derived_params.size(); ++i ){
        datalist.ParList->Add_Derived_Param(derived_params[i]);
        // cout << "adding " << derived_params[i] << endl;
    }

    // datalist.PrintParamList();
    // exit(0);
}


bool CosmoTheory::update( imcmc_double& full_param ) {
    bool status = engine->updateParValues(class_mcmc_params, full_param);
    return status;
}

void CosmoTheory::add_derived_params( string& paramfile ) {

    has_derived_params = false; // by default assumes no derived parameters.

    if( Read::Has_Key_in_File(paramfile, "derived_cosmo_params") ) {

        int nderived = Read::Num_of_Value_for_Key(paramfile, "derived_cosmo_params");

        if( nderived > 0 ) {

            string* pnames = new string[nderived];
            Read::Read_Array_of_String_from_File( paramfile, "derived_cosmo_params", pnames, nderived );

            for( int i=0; i<nderived; ++i ) {
                derived_params.push_back(pnames[i]);
// std::cout << " add derived cosmo param: " << pnames[i] << "\n";
            }

            delete[] pnames;

            has_derived_params = true;
        }
    }
}

void CosmoTheory::update_derived_params( imcmc_double& full_param ) {

    //  this function can be further improved, remove the while loop.

    imcmc_vector_string_iterator it = derived_params.begin();

    while( it != derived_params.end() ) {

		//cout << "updating derived param: " << *it << endl;

        if( *it == "100*theta_s" ) {
            full_param[*it] = engine->get_100theta_s();
        }

        if( *it == "theta_s" ) {
            full_param[*it] = engine->get_100theta_s()/100.0;
        }
        
        if( *it == "tau_reio" ) {
            full_param[*it] = engine->getTauReio();
        }
        
        if( *it == "H0" ) {
            full_param[*it] = engine->get_H0();
        }
        
        if( *it == "Omega_m" ) {
            full_param[*it] = engine->get_Omega_m();
        }
        
        if( *it == "omegamh2") {
            double H0 = engine->get_H0();
            double Omega_m = engine->get_Omega_m();
            full_param[*it] = Omega_m*(H0/100.0)*(H0/100.0);
        }

        if( *it == "Omega_b" ) {
            full_param[*it] = engine->get_Omega_b();
        }

        if( *it == "omegabh2") {
            double H0 = engine->get_H0();
            double Omega_b = engine->get_Omega_b();
            full_param[*it] = Omega_b*(H0/100.0)*(H0/100.0);
        }
        
        if( *it == "lnAs") {

            double A_s = -1.0;
            if( opt_As == 0 )
                A_s = full_param["A_s"];
            else if( opt_As == 1 )
                A_s = full_param["10^9A_s"]*1e-9;
            else if( opt_As == 2 )
                A_s = exp(full_param["ln(10^10A_s)"])*1e-10;

            if( A_s < 0 ) {
                MPI_ClassMC_ERROR("failed to get A_s!");
            }
            full_param[*it] = log(A_s);
        }
        
        if( (*it == "YHe_bbn") && do_CMB ) {
            if( engine->YHe_bbn > 0 )
                full_param[*it] = engine->YHe_bbn;
            else {//  This should never happen!!!
                MPI_ClassMC_ERROR("#  failed to get engine->YHe_bbn.");
            }
        }

        if( *it == "sigma8" ){
            full_param[*it] = engine->get_sigma8(0.0);  // get sigma8 at z=0
//			cout << "sigma8 = " << full_param[*it] << endl;
        }
        
        // =============================================================
        // add other derieved parameters here ...
        //

        if( *it == "DDE_w_z2.5" ){
        	double w, weff;
        	engine->get_EoS(2.4995, &w, &weff);
            full_param[*it] = w;
        }

        ++it;
    }
}
