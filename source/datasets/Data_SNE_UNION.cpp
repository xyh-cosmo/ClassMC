/*
	TODO: Update data reading part, can use armadillo's member functions.
*/

#include "SNE.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
using namespace arma;

void SNe_UNION::Init( string& dataset ){

    data_info.GetInfo(dataset);

	std::string data 	= Read::Read_String_from_File( dataset, "data_file" );
	int size 			= Read::Read_Int_from_File( dataset, "data_size" );
	sne_num = size;

	z   = arma::zeros(size, 1);
	mu  = arma::zeros(size, 1);
	dmu = arma::zeros(size, 1);
	P   = arma::zeros(size, 1);
	icov= arma::zeros(size, size);

	systematic = Read::Read_Bool_from_File( dataset, "include_systematics" );
	string covmat;
	if( systematic == true ){
		cout << "==> union2.1 using the covariance matrix with systematics\n";
		covmat = Read::Read_String_from_File( dataset, "cov_sys" );
	}
	else if( systematic == false ){
		cout << "==> union2.1 using the covariance matrix without systematics\n";
		covmat = Read::Read_String_from_File( dataset, "cov_nosys" );
	}

//	now read in data
	std::ifstream infile_data(data.c_str());
	std::string line, sn_name;

	int i=0;
	while( std::getline(infile_data, line) ){
	    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
			std::stringstream stream(line);
			stream >> sn_name >> z(i) >> mu(i) >> dmu(i) >> P(i);
			++i;
	    }
	}

	infile_data.close();

	if( i != size )
	    throw runtime_error("*** UNION::ReadData() ==> wrong number of UNION2.1 distance moduli data !");

//  now let's read in covariance matrix

	std::ifstream infile_cov(covmat.c_str());	// size is 580*580

	if( !infile_cov )
	    throw runtime_error("*** UNION::ReadData() ==> cannot open: " + covmat);
	else{
		i = 0;
		while( std::getline(infile_cov, line) ){
	        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
				std::stringstream stream(line);
				double C_ij;
				for(int j=0; j<size; ++j){
					stream >> C_ij;
	                icov(i,j) = C_ij;
				}
				++i;
	        }
		}
	}

	infile_cov.close();

	if( i != size )
	    throw runtime_error("*** UNION::ReadData() ==> wrong number of UNION2.1 distance moduli cov data !");
}
