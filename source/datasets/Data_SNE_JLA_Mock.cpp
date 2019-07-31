//  ========================================================
//  ==================   SNLS   ============================
//  ========================================================

#include "SNE.hpp"
#include "Misc.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

SNe_JLA_Mock::SNe_JLA_Mock(){
	sne_num = 0;
	has_err = true;
	use_full_cov = true;
	invert_err = false;
}

SNe_JLA_Mock::~SNe_JLA_Mock(){
	// do nothing!
}

void SNe_JLA_Mock::Init( std::string& dataset ){
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

	cout << "==> Loading SNe_JLA_Mock data from: " << data_file << endl;

	z 	= arma::vec(size,1);
	mb 	= arma::vec(size,1);
	dmb = arma::vec(size,1);
	arma::vec mb0 = arma::vec(size,1);

	string line, sn_name;
	ifstream infile(data_file.c_str());

	if( !infile.good() ){
		MPI_ClassMC_ERROR("failed to open: "+data_file);
	}

	int i = 0;
	double temp;
	while( getline(infile,line) ){
	// TODO: remove SNe names (strings)
		if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
			stringstream stream(line);
			if( has_err )
				stream >> sn_name >> z(i) >> mb(i) >> dmb(i) >> mb0(i);
			else
				stream >> sn_name >> z(i) >> temp >> dmb(i) >> mb(i);
			++i;
			if( i > size ){
				MPI_ClassMC_ERROR("# error in loading JLA mock!\n");
			}
		}
	}

	infile.close();

	cout << "==> Loading SNe_JLA_Mock inverse covariance matrix from: " << icov_file << endl;
	covmat_inv.load(icov_file,arma::raw_ascii);
	if( int(covmat_inv.n_rows) != size ){
		MPI_ClassMC_ERROR("# error in loading JLA mock (inv-) covariance matrix!\n");
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


