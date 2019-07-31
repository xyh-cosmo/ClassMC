#include "SNE.hpp"
#include "Misc.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
// using namespace arma;

SNe_WFIRST::SNe_WFIRST(){
	cout << "==> Creating SNe_WFIRST\n";
}

SNe_WFIRST::~SNe_WFIRST(){
	cout << "==> Destroying SNe_WFIRST\n";
}

void SNe_WFIRST::Init( string& dataset ){

	int size;
	string data;

	data_info.GetInfo(dataset);

	data = Read::Read_String_from_File(dataset, "data_file");
	size = Read::Read_Int_from_File(dataset, "size");

	has_err = true;
	if( Read::Has_Key_in_File(dataset,"has_err") )
		has_err = Read::Read_Bool_from_File(dataset,"has_err");

	invert_err = false;
	if( Read::Has_Key_in_File(dataset,"invert_err") )
		invert_err = Read::Read_Bool_from_File(dataset,"invert_err");

	arma::mat wfirst;

	wfirst.load(data,arma::raw_ascii);  // read in as a matrix

	if( size != int(wfirst.n_rows) ){
		string errmsg = "# failed to load WFIRST mock data, size = " + Read::IntToString(size);
		errmsg += " does not match the number of rows of : " + data + "\n";
		MPI_ClassMC_ERROR(errmsg);		
//		throw runtime_error(errmsg);
	}
	else{
		sne_num = size;
		z	= arma::zeros(size,1);
		mu	= arma::zeros(size,1);
		dmu	= arma::zeros(size,1);
		arma::vec mu0 = arma::vec(size,1);

		if( has_err ){
			for( int i=0; i<size; ++i ){
				z(i)	= wfirst(i,0);
				mu(i)	= wfirst(i,1);
				dmu(i)	= wfirst(i,2);
				mu0(i)  = wfirst(i,3);
			}
		}
		else{
			for( int i=0; i<size; ++i ){
				z(i)	= wfirst(i,0);
				mu(i)	= wfirst(i,3);
				dmu(i)	= wfirst(i,2);
			}
		}

		if( invert_err && has_err ){
			mu = mu0 - (mu-mu0);
		}
	}
}
