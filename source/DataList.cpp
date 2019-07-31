#include "DataList.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

DataList::DataList() {

    ParList = NULL;

//  data option flags
    use_Age = false;
    use_BAO = false;
    use_HST	= false;
    use_Hz	= false;
    use_RSD = false;

    use_CMB		= false;
#if defined(_USE_WMAP7_)
    use_WMAP7   = false;
#elif defined(_USE_WMAP9_)
	use_WMAP9	= false;
#endif
    use_PLK		= false;

    use_SN			= false;
    use_UNION		= false;
    use_SNLS		= false;
    use_SNLS_Mock	= false;
    use_JLA			= false;
    use_JLA_Mock	= false;
	use_WFIRST		= false;

    use_DDE_CPZ_prior	= false;
    use_CMB_DistPrior   = false;

//  data structure pointers
    data_age 	= NULL;
    data_bao 	= NULL;
    data_HST    = NULL;
    data_Hz     = NULL;

    data_rsd    = NULL;

    data_sne_union      = NULL;
    data_sne_snls       = NULL;
    data_sne_snls_mock  = NULL;
    data_sne_jla        = NULL;
    data_sne_jla_mock   = NULL;
	data_sne_wfirst		= NULL;

    data_plk2015        = NULL;

//  prior structure pointer
    prior_dde_cpz   = NULL;
    prior_cmb_dist  = NULL;
}

DataList::~DataList() {

    if( ParList != NULL ) {
        ParList->~ParamList();
        ParList = NULL;
    }

    if( use_Age && data_age != NULL ){
    	data_age->~Data_Age();
    	data_age = NULL;
    }

    if( use_BAO && data_bao != NULL ){
    	data_bao->~Data_BAO();
    	data_bao = NULL;
    }

    if( use_HST && data_HST != NULL ) {
        data_HST->~Data_HST();
        data_HST = NULL;
    }

    if( use_Hz && data_Hz != NULL ) {
        data_Hz->~Data_Hubble();
        data_Hz = NULL;
    }

    if( use_RSD && data_rsd != NULL ){
        data_rsd->~Data_RSD();
        data_rsd = NULL;
    }

    if( use_CMB ) {
        if( use_PLK && data_plk2015 != NULL ) {
            data_plk2015->~Data_Planck2015();
            data_plk2015 = NULL;
        }
#if defined(_USE_WMAP7_)
		else if( use_WMAP7 && data_wmap7 != NULL ){
			data_wmap7->~Data_WMAP7();
			data_wmap7 = NULL;
		}
#elif defined(_USE_WMAP9_)
        else if( use_WMAP9 && data_wmap9 != NULL ){
			data_wmap9->~Data_WMAP9();
			data_wmap9 = NULL;
		}
#endif
    }

    if( use_SN ) {
        if( use_UNION && data_sne_union != NULL ) {
            data_sne_union->~SNe_UNION();
            data_sne_union = NULL;
        }

        if( use_SNLS && data_sne_snls != NULL ) {
            data_sne_snls->~SNe_SNLS();
            data_sne_snls = NULL;
        }

        if( use_SNLS_Mock && data_sne_snls_mock != NULL ){
            data_sne_snls_mock->~SNe_SNLS_Mock();
            data_sne_snls_mock = NULL;
        }

        if( use_JLA && data_sne_jla != NULL ) {
            data_sne_jla->~SNe_JLA();
            data_sne_jla = NULL;
        }

        if( use_JLA_Mock && data_sne_jla_mock != NULL ){
            data_sne_jla_mock->~SNe_JLA_Mock();
            data_sne_jla_mock = NULL;
        }

		if( use_WFIRST && data_sne_wfirst != NULL ) {
			data_sne_wfirst->~SNe_WFIRST();
			data_sne_wfirst = NULL;
		}

    }

    if( use_DDE_CPZ_prior && prior_dde_cpz != NULL ){
        prior_dde_cpz->~DDE_CPZ();
        prior_dde_cpz = NULL;
    }

    if( use_CMB_DistPrior && prior_cmb_dist != NULL ){
        prior_cmb_dist->~CMB_Dist();
        prior_cmb_dist = NULL;
    }
}

void DataList::Init( string& paramfile ) {

    MPI_cout("Start to configure DataList from input file : " + paramfile );

    ParList = new ParamList;

    vector<string> used_data;

    use_Age = Read::Read_Bool_from_File(paramfile, "use_Age");
    if( use_Age ){

        string Age_dataset = "";
        Age_dataset = Read::Read_String_from_File(paramfile,"Age_dataset");

        if( Age_dataset == "" ){
            MPI_ClassMC_ERROR("==> failed to read Age_dataset");
        }

        data_age = new Data_Age;
        data_age->Init(Age_dataset);
        ParList->Add_Params(data_age->data_info);
        used_data.push_back("Galaxy Age");
    }

    use_BAO = Read::Read_Bool_from_File(paramfile, "use_BAO");
    if( use_BAO ) {

        string BAO_dataset = "";
        BAO_dataset = Read::Read_String_from_File(paramfile,"BAO_dataset");

        if( BAO_dataset == "" ) {
            MPI_ClassMC_ERROR("==> failed to get BAO_dataset");
        }
        
		data_bao = new Data_BAO;
        data_bao->Init(BAO_dataset);
        ParList->Add_Params(data_bao->data_info);
        used_data.push_back(data_bao->data_info.DataName);
    }

    use_HST = Read::Read_Bool_from_File(paramfile, "use_HST");
    if( use_HST ) {

        string HST_dataset = "";
        HST_dataset = Read::Read_String_from_File(paramfile,"HST_dataset");

		if( HST_dataset == "" ){
			MPI_ClassMC_ERROR("==> failed to get HST_dataset");
		}

        data_HST = new Data_HST;
        data_HST->Init(HST_dataset);
        ParList->Add_Params(data_HST->data_info);
        used_data.push_back("HST H0");
    }

    use_Hz = Read::Read_Bool_from_File(paramfile, "use_Hz");
    if( use_Hz ) {

        string Hz_dataset = "";
        Hz_dataset = Read::Read_String_from_File(paramfile,"Hz_dataset");
        
		if( Hz_dataset == "" ){
            MPI_ClassMC_ERROR("==> failed to get Hz_dataset");
        }

        data_Hz = new Data_Hubble;
        data_Hz->Init(Hz_dataset);
        ParList->Add_Params(data_Hz->data_info);
        used_data.push_back("Hubble parameter");
    }

    use_RSD = Read::Read_Bool_from_File(paramfile, "use_RSD");
    if( use_RSD ){
        string rsd_dataset = "";
        rsd_dataset = Read::Read_String_from_File(paramfile,"RSD_dataset");

        if( rsd_dataset == "" ){
            MPI_ClassMC_ERROR("===> failed to get rsd_dataset from :"+paramfile);
        }

        data_rsd = new Data_RSD;
        data_rsd->Init(rsd_dataset);
        ParList->Add_Params(data_rsd->data_info);
        used_data.push_back(data_rsd->data_info.DataName);
    }

    use_CMB = Read::Read_Bool_from_File(paramfile, "use_CMB");
    if( use_CMB == true ) {

#if defined(_USE_WMAP7_)
        use_WMAP7 = Read::Read_Bool_from_File(paramfile, "use_WMAP7");
#elif defined(_USE_WMAP9_)
		use_WMAP9 = Read::Read_Bool_from_File(paramfile, "use_WMAP9");
#endif
        use_PLK = Read::Read_Bool_from_File(paramfile, "use_PLK");

#if defined(_USE_WMAP7_)
		if( use_WMAP7 && use_PLK ){
			MPI_ClassMC_ERROR("==> Error, only one of WMAP7 and Planck can be used at one time");
		}
#elif defined(_USE_WMAP9_)
		if( use_WMAP9 && use_PLK ){
			MPI_ClassMC_ERROR("==> Error, only one of WMAP9 and Planck can be used at one time");
		}
#endif

        if( use_PLK ) {
            string plk_dataset = "";
            plk_dataset = Read::Read_String_from_File(paramfile,"plk_dataset");

            if( plk_dataset=="" ){
                MPI_ClassMC_ERROR("==> failed to get plk_dataset");
            }

            data_plk2015 = new Data_Planck2015;
            data_plk2015->Init(plk_dataset);
            ParList->Add_Params(data_plk2015->data_info);
            used_data.push_back("Planck CMB");
        }
#if defined(_USE_WMAP7_)
		else if( use_WMAP7 ){
			string wmap7_dataset = "";	// not used acutally ...
			wmap7_dataset = Read::Read_String_from_File(paramfile,"wmap7_dataset");

			if( wmap7_dataset == "" ){
				MPI_ClassMC_ERROR("==> failed to get wmap7_dataset");
			}

			data_wmap7 = new Data_WMAP7;
			data_wmap7->Init(wmap7_dataset);
			ParList->Add_Params(data_wmap7->data_info);
            used_data.push_back("WAMP7 CMB");
		}
#elif defined(_USE_WMAP9_)
        else if( use_WMAP9 ){
			string wmap9_dataset = "";	// not used acutally ...
			wmap9_dataset = Read::Read_String_from_File(paramfile,"wmap9_dataset");

			if( wmap9_dataset == "" ){
				MPI_ClassMC_ERROR("==> failed to get wmap9_dataset");
			}

			data_wmap9 = new Data_WMAP9;
			data_wmap9->Init(wmap9_dataset);
			ParList->Add_Params(data_wmap9->data_info);
            used_data.push_back("WAMP9 CMB");
		}
#endif
    }

    use_SN = Read::Read_Bool_from_File(paramfile, "use_SN");
    if( use_SN ) {

        use_UNION       = Read::Read_Bool_from_File(paramfile, "use_UNION");
        use_SNLS        = Read::Read_Bool_from_File(paramfile, "use_SNLS");
        use_SNLS_Mock   = Read::Read_Bool_from_File(paramfile, "use_SNLS_Mock");
        use_JLA         = Read::Read_Bool_from_File(paramfile, "use_JLA");
        use_JLA_Mock    = Read::Read_Bool_from_File(paramfile, "use_JLA_Mock");
		use_WFIRST		= Read::Read_Bool_from_File(paramfile, "use_WFIRST");

        int cnt = 0;
        if( use_UNION == true ) ++cnt;
        if( use_SNLS == true ) ++cnt;
        if( use_SNLS_Mock == true ) ++cnt;
        if( use_JLA == true ) ++cnt;
        if( use_JLA_Mock == true ) ++cnt;

		// note WFIRST can be used with any of the above SNe sample (real or mock)

        if( cnt > 1 ){
            MPI_ClassMC_ERROR("Only one of {Union, SNLS3, SNLS3_Mock, JLA, JLA_Mock} can be used at one time!");
        }

        if( use_UNION ) {

            string union_dataset = "";
            union_dataset = Read::Read_String_from_File(paramfile,"UNION_dataset");

            if( union_dataset=="" ){
				MPI_ClassMC_ERROR("==> failed to get union_dataset");
            }

            data_sne_union = new SNe_UNION;
            data_sne_union->Init(union_dataset);
            ParList->Add_Params(data_sne_union->data_info);
            used_data.push_back("Union2.1");
        }

        if( use_SNLS ) {

            string snls_dataset = "";
            snls_dataset = Read::Read_String_from_File(paramfile,"SNLS_dataset");

            if( snls_dataset=="" ){
				MPI_ClassMC_ERROR("failed to get snls_dataset");
            }

            data_sne_snls = new SNe_SNLS;
            data_sne_snls->Init(snls_dataset);
            ParList->Add_Params(data_sne_snls->data_info);
            used_data.push_back("SNLS3");
        }

        if( use_SNLS_Mock ){
            string snls_dataset = "";
            snls_dataset = Read::Read_String_from_File(paramfile,"SNLS_Mock_dataset");
            
            if( snls_dataset=="" ){
                MPI_ClassMC_ERROR("==> failed to get \'snls_dataset\'");
            }

            data_sne_snls_mock = new SNe_SNLS_Mock;
            data_sne_snls_mock->Init(snls_dataset);
            ParList->Add_Params(data_sne_snls_mock->data_info);
            used_data.push_back("SNLS3_Mock");
        }

        if( use_JLA ) {

            string jla_dataset = "";
            jla_dataset = Read::Read_String_from_File(paramfile,"JLA_dataset");

            if( jla_dataset=="" ){
				MPI_ClassMC_ERROR("failed to get jla_dataset");
            }

            data_sne_jla = new SNe_JLA(0);
            data_sne_jla->Init(jla_dataset);
            ParList->Add_Params(data_sne_jla->data_info);
            used_data.push_back("JLA");
            
            //	debug
//            cout << "successfully loaded JLA SNe Ia data ...\n";
        }

        if( use_JLA_Mock ){
			string jla_mock_dataset = "";
			jla_mock_dataset = Read::Read_String_from_File(paramfile,"JLA_Mock_dataset");
			
			if( jla_mock_dataset == "" ){
				MPI_ClassMC_ERROR("==> failed to get jla_Mock_dataset");
			}
			
			data_sne_jla_mock = new SNe_JLA_Mock;
			data_sne_jla_mock->Init(jla_mock_dataset);
			ParList->Add_Params(data_sne_jla_mock->data_info);
			used_data.push_back("JLA_Mock");
        }

		if( use_WFIRST ){
			string wfirst_dataset = "";
			wfirst_dataset = Read::Read_String_from_File(paramfile,"wfirst_dataset");
			if( wfirst_dataset == "" ){
				MPI_ClassMC_ERROR("==> failed to get wfirst_dataset");
			}

			data_sne_wfirst = new SNe_WFIRST;
			data_sne_wfirst->Init(wfirst_dataset);
			ParList->Add_Params(data_sne_wfirst->data_info);
			used_data.push_back("WFIRST_Mock");
		}

    }

    // for debug
    // PrintParamList();
    // exit(0);

    //  prior can be viewed as some kind of dataset
    use_DDE_CPZ_prior = Read::Read_Bool_from_File(paramfile,"use_DDE_CPZ_prior");
    if( use_DDE_CPZ_prior ){

        string DDE_CPZ_prior_settings = "";
        DDE_CPZ_prior_settings = Read::Read_String_from_File(paramfile,"DDE_CPZ_prior_settings");

        if( DDE_CPZ_prior_settings=="" ){
            MPI_ClassMC_ERROR("==> failed to get DDE_CPZ_prior_settings");
        }

        prior_dde_cpz = new DDE_CPZ;
        prior_dde_cpz->Init(DDE_CPZ_prior_settings);
        ParList->Add_Params(prior_dde_cpz->data_info);
        used_data.push_back("DDE CPZ prior");
    }

    use_CMB_DistPrior = Read::Read_Bool_from_File(paramfile,"use_CMB_DistPrior");
    if( use_CMB_DistPrior ){

        if( use_CMB ){
			MPI_ClassMC_ERROR("==> full CMB data and CMB distance prior cannot be used together!");
        }

        string CMB_distance_prior_dataset = "";
        CMB_distance_prior_dataset = Read::Read_String_from_File(paramfile,"CMB_distance_prior_dataset");

        if( CMB_distance_prior_dataset=="" ){
            MPI_ClassMC_ERROR("==> failed to get CMB_distance_prior_dataset");
        }

        prior_cmb_dist = new CMB_Dist;
        prior_cmb_dist->Init(CMB_distance_prior_dataset);
        ParList->Add_Params(prior_cmb_dist->data_info);
        used_data.push_back("CMB distance prior: "+prior_cmb_dist->data_info.DataName);
    }

	if( MPI::COMM_WORLD.Get_rank() == ROOT_RANK ){
		cout << "==> Finished DataList configuration.\n";
        cout << "==> used datasets and prior:\n";

        string chain_root;
        if( Read::Has_Key_in_File(paramfile, "chain_root") )
            chain_root = Read::Read_String_from_File(paramfile, "chain_root");
        else
            chain_root = paramfile;
        
        string fname = chain_root + ".data_and_prior";
        ofstream save_used_data(fname.c_str());
        save_used_data << "# used data and prior\n";

        for( size_t i=0; i<used_data.size(); ++i ){
            cout << " -- " << used_data[i] << "\n";
            save_used_data << "dataset_" << i+1 << " = " << used_data[i] << "\n";
        }

        save_used_data.close();
	}
	
	// decide whether to compute sigma8
	compute_sigma8 = false;
	if( Read::Has_Key_in_File(paramfile, "compute_sigma8") ){
	    compute_sigma8 = Read::Read_Bool_from_File(paramfile, "compute_sigma8");
	}
}

void DataList::PrintParamList() {
    cout << "\n##### Caution:\n This function needs improments #####\n";
    ParList->Print_MCMC_Params();
    ParList->Print_Derived_Params();
}

void DataList::Get_MCMC_Params( vector<string>& mcmc_params ) {

    mcmc_params.clear();
    map<string,string>::iterator it;

    for( it=ParList->MCMC_Params.begin(); it!=ParList->MCMC_Params.end(); ++it ) {
        mcmc_params.push_back(it->first);
    }
}

void DataList::Get_Derived_Params( vector<string>& derived_params ) {

    derived_params.clear();
    map<string,string>::iterator it;

    for( it=ParList->Derived_Params.begin(); it!=ParList->Derived_Params.end(); ++it ) {
        derived_params.push_back(it->first);
    }
}
