//  ========================================================
//  ==================   SNLS   ============================
//  ========================================================

#include "SNE.hpp"
#include "Misc.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

SNe_SNLS::SNe_SNLS(){
	// do nothing
    marg_scriptm = false;
}

SNe_SNLS::~SNe_SNLS(){
	// do nothing
}

void SNe_SNLS::ReadCov( arma::mat& covmat, string& covmat_file ){
    double covx;
    string line;

    ifstream file(covmat_file.c_str());

    if( !file ) throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + covmat_file);

    int I = 0, J = 0;
    getline(file, line);  // first line is just the number 472
    while( getline(file, line) ){
        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
            std::stringstream stream(line);
            while( stream >> covx ){
                if( J < sne_num ){
                    covmat(I,J) = covx;
                    ++J;
                }
            }
            if( J == sne_num ){
                ++I;
                J = 0;
            }
        }
    }

    file.close();
}

void SNe_SNLS::Init( string& dataset ){

    data_info.GetInfo(dataset);
	string data_file    = Read::Read_String_from_File(dataset, "data_file");
	int size            = Read::Read_Int_from_File(dataset, "size");
	sne_num = size;

//  allocate memory
//	1)	light-curve data
	zcmb        = arma::zeros(size,1);
	zhel        = arma::zeros(size,1);
	dz          = arma::zeros(size,1);
	mb          = arma::zeros(size,1);
	dmb         = arma::zeros(size,1);
	s           = arma::zeros(size,1);
	ds          = arma::zeros(size,1);
	c           = arma::zeros(size,1);
	dc          = arma::zeros(size,1);
	var3        = arma::zeros(size,1);
	dvar3       = arma::zeros(size,1);
	cov_m_s     = arma::zeros(size,1);
	cov_m_c     = arma::zeros(size,1);
	cov_s_c     = arma::zeros(size,1);

//  pre-computed variables
    pre_vars    = arma::zeros(size,1);

//	2)	covariance matrix
	cov_mB_mB       = arma::zeros(size,size);
	cov_mB_alpha    = arma::zeros(size,size);
	cov_mB_beta     = arma::zeros(size,size);
	cov_alpha_alpha = arma::zeros(size,size);
	cov_alpha_beta  = arma::zeros(size,size);
	cov_beta_beta   = arma::zeros(size,size);

	cov_tot         = arma::zeros(size,size);
	icov_tot        = arma::zeros(size,size);

	pecz            = Read::Read_Double_from_File(dataset, "pecz");

	intrinsicdisp   = Read::Read_Double_from_File(dataset, "intrinsicdisp");

	use_four_disp	= Read::Read_Bool_from_File(dataset, "use_four_disp");
	intrinsicdisp0  = Read::Read_Double_from_File(dataset, "intrinsicdisp0");
	intrinsicdisp1  = Read::Read_Double_from_File(dataset, "intrinsicdisp1");
	intrinsicdisp2  = Read::Read_Double_from_File(dataset, "intrinsicdisp2");
	intrinsicdisp3  = Read::Read_Double_from_File(dataset, "intrinsicdisp3");

	twoscriptmfit   = Read::Read_Bool_from_File(dataset, "twoscriptmfit");
    marg_scriptm    = Read::Read_Bool_from_File(dataset, "marg_scriptm");
	scriptmcut      = Read::Read_Double_from_File(dataset, "scriptmcut");

	string mag_covmat_file;
	string mag_stretch_covmat_file;
	string mag_colour_covmat_file;
	string stretch_covmat_file;
	string stretch_colour_covmat_file;
	string colour_covmat_file;

	mag_covmat_file             = Read::Read_String_from_File(dataset, "mag_covmat_file");
	stretch_covmat_file         = Read::Read_String_from_File(dataset, "stretch_covmat_file");
	colour_covmat_file          = Read::Read_String_from_File(dataset, "colour_covmat_file");
	mag_stretch_covmat_file     = Read::Read_String_from_File(dataset, "mag_stretch_covmat_file");
	mag_colour_covmat_file      = Read::Read_String_from_File(dataset, "mag_colour_covmat_file");
	stretch_colour_covmat_file  = Read::Read_String_from_File(dataset, "stretch_colour_covmat_file");

//	read data
	std::string line;
	int sn_num=0;

	std::ifstream infile( data_file.c_str() );

	if( !infile )
	    throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + data_file);

	std::getline(infile, line);    // skip comments in the first line
	while( std::getline(infile, line) ){
	    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		    std::string sn_namex;
		    double zcmbx;
		    double zhelx;
		    double dzx;
		    double mbx;
		    double dmbx;
		    double sx;
		    double dsx;
		    double cx;
		    double dcx;
		    double var3x;
		    double dvar3x;
		    double cov_m_sx;
		    double cov_m_cx;
		    double cov_s_cx;
		    int setx;

			std::stringstream stream(line);
			stream 	>> sn_namex
					>> zcmbx >> zhelx >> dzx
					>> mbx >> dmbx
					>> sx >> dsx
					>> cx >> dcx
					>> var3x >> dvar3x
					>> cov_m_sx
					>> cov_m_cx
					>> cov_s_cx
					>> setx;

			sne_name.push_back(sn_namex);
			set.push_back(setx);

	        zcmb(sn_num)    = zcmbx;
	        zhel(sn_num)    = zhelx;
	        dz(sn_num)      = dzx;
	        mb(sn_num)      = mbx;
	        dmb(sn_num)     = dmbx;
	        s(sn_num)       = sx;
	        ds(sn_num)      = dsx;
	        c(sn_num)       = cx;
	        dc(sn_num)      = dcx;
	        var3(sn_num)    = var3x;
	        dvar3(sn_num)   = dvar3x;
	        cov_m_s(sn_num) = cov_m_sx;
	        cov_m_c(sn_num) = cov_m_cx;
	        cov_s_c(sn_num) = cov_s_cx;

			++sn_num;
	    }
	}

	infile.close();

	if( sn_num != size )
	    throw runtime_error("*** SNLS::ReadData() ==> Fatal error: wrong number of SNe light-curve data!");

    if( marg_scriptm && twoscriptmfit ){

        K1 = arma::zeros(size,1);
        K2 = arma::zeros(size,1);

        for( int i=0; i<sne_num; ++i ){
			if( var3(i) > scriptmcut ){
            	K1(i) = 0;
            	K2(i) = 1;
			}
			else{
				K1(i) = 1;
				K2(i) = 0;
			}
        }
    }

    ReadCov( cov_mB_mB, mag_covmat_file );

//  stretch_covmat
    ReadCov( cov_alpha_alpha, stretch_covmat_file );

//  colour_covmat
    ReadCov( cov_beta_beta, colour_covmat_file );

//  mag_stretch
    ReadCov( cov_mB_alpha, mag_stretch_covmat_file );

//  mag_colour
    ReadCov( cov_mB_beta, mag_colour_covmat_file );

//  stretch_colour
    ReadCov( cov_alpha_beta, stretch_colour_covmat_file );


//  fill pre_vars
    double zfacsq   = ( 5.0/log(10.0) )*( 5.0/log(10.0) );
    double intdisp[4] = {intrinsicdisp0,intrinsicdisp1,intrinsicdisp2,intrinsicdisp3};
    double emptyfac;

    for(int i=0; i<sne_num; ++i){
    //  light curve fitting uncertainties
        pre_vars(i) = dmb(i) * dmb(i);
    //  peculiar velocity and redshift error...
        emptyfac = ( 1.0 + zcmb(i) ) / ( zcmb(i) * (1 + 0.5*zcmb(i)));
        pre_vars(i) += ( dz(i)*dz(i) + pecz*pecz )*zfacsq*emptyfac*emptyfac;
    //  intrinsic dispersion
        if( use_four_disp == true )
            pre_vars(i) += (intdisp[set[i]]*intdisp[set[i]]);
        else
            pre_vars(i) += (intrinsicdisp * intrinsicdisp);
        pre_vars(i) += pow(0.055*zcmb(i), 2);   //  lensing: sigma_lens = 0.055*z
        // CosmoMC didnot include this lensing term, WHY?
    }

}

void SNe_SNLS::UpdateNuisance( imcmc_double par ){
	alpha   = par["alpha_snls"];
	beta    = par["beta_snls"];

    if( !marg_scriptm ){
        if( !twoscriptmfit )
            scriptMB = par["scriptMB_snls"];
        else{
            scriptMB1 = par["scriptMB1_snls"];
            scriptMB2 = par["scriptMB2_snls"];
        }
    }
}

void SNe_SNLS::UpdateCov(){

	// cov_tot.fill(0);
    double alpha2 = alpha*alpha;
    double beta2 = beta*beta;
    double alphabeta = alpha*beta;

	cov_tot = ( cov_mB_mB
	        +   alpha2*cov_alpha_alpha
	        +   beta2*cov_beta_beta
	        +   2*alpha*cov_mB_alpha
	        -   2*beta*cov_mB_beta
	        -   2*alphabeta*cov_alpha_beta );

//  add diagonal terms
	for(int i=0; i<sne_num; ++i){
	    cov_tot(i,i)   += ( pre_vars(i)
	                    +   alpha2 * ds(i) * ds(i)
	                    +   beta2 * dc(i) * dc(i)
	                    +   2*alpha * cov_m_s(i)
	                    -   2*beta * cov_m_c(i)
	                    -   2*alphabeta * cov_s_c(i) );
	}

//	now invert the covariance matrix
	icov_tot = cov_tot.i();
}

void SNe_SNLS::SaveTotalCov( std::string tot_cov_rootname ){
	string cov_fname = tot_cov_rootname + ".txt";
	string icov_fname= tot_cov_rootname + "_inv.txt";
	cov_tot.save(cov_fname,arma::raw_ascii);
	icov_tot.save(icov_fname,arma::raw_ascii);
}


SNe_SNLS_Mock::SNe_SNLS_Mock(){
	sne_num = 0;
	has_err = false;
	use_full_cov = true;
	invert_err = false;
}

SNe_SNLS_Mock::~SNe_SNLS_Mock(){
	// do nothing
}


void SNe_SNLS_Mock::Init( std::string& dataset ){
	int size;
	string data_file, icov_file;
	data_info.GetInfo(dataset);

	if( Read::Has_Key_in_File(dataset,"data_file"))
		data_file = Read::Read_String_from_File(dataset,"data_file");
	else
		MPI_ClassMC_ERROR("cannot find \'data_file\'");

	if( Read::Has_Key_in_File(dataset,"has_err") )
		has_err = Read::Read_Bool_from_File(dataset,"has_err");

	if( Read::Has_Key_in_File(dataset,"invert_err") )
		invert_err = Read::Read_Bool_from_File(dataset,"invert_err");

	if( (has_err==false) && (invert_err==true) ){
		MPI_ClassMC_ERROR("invert_err can be used only when has_err==true");
	}

	if( Read::Has_Key_in_File(dataset,"use_full_cov") )
		use_full_cov = Read::Read_Bool_from_File(dataset,"use_full_cov");

	icov_file = Read::Read_String_from_File(dataset,"icov_file");
	size = Read::Read_Int_from_File(dataset,"size");

	cout << "==> Loading SNe_SNLS_Mock data from: " << data_file << endl;

	z 	= arma::vec(size,1);
	mb 	= arma::vec(size,1);
	dmb = arma::vec(size,1);
	arma::vec mb0 = arma::vec(size,1);

	string line, sn_name;
	ifstream infile(data_file.c_str());

	int i = 0;
	double temp;
	while( getline(infile,line) ){
	// TODO: remove the SNe names (strings)
		if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
			stringstream stream(line);
			if( has_err )
				stream >> sn_name >> z(i) >> mb(i) >> dmb(i) >> mb0(i);
			else
				stream >> sn_name >> z(i) >> temp >> dmb(i) >> mb(i);
			++i;
			if( i > size ){
				MPI_ClassMC_ERROR("# error in loading SNLS3 mock!\n");
			}
		}
	}

	infile.close();

	cout << "==> Loading SNe_SNLS_Mock inverse covariance matrix from: " << icov_file << endl;
	covmat_inv.load(icov_file,arma::raw_ascii);
	if( int(covmat_inv.n_rows) != size ){
		MPI_ClassMC_ERROR("# error in loading SNLS3 mock (inv-) covariance matrix!\n");
	}
	// exit(0);

	if( use_full_cov == false ){
		arma::mat covmat_inv_temp = covmat_inv;
		covmat_inv.fill(0.0);
		covmat_inv.diag() = covmat_inv_temp.diag();
	}

	sne_num = size;

	if( invert_err == true ){
		mb = mb0 - (mb-mb0);
	}
}


