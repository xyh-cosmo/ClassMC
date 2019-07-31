#include <imcmc/parser++.hpp>
#include "Age.hpp"

using namespace std;
using namespace imcmc::parser;


Data_Age::Data_Age(){
	age_num	= 0;
	z		= NULL;
	age		= NULL;
	age_err	= NULL;
	age_init= false;
	gft_marg= false;
}

Data_Age::~Data_Age(){
	if( (age_num > 0) && age_init ){
		delete[] z;
		delete[] age;
		delete[] age_err;
	}
}

void Data_Age::Init( std::string& age_dataset ){

	data_info.GetInfo(age_dataset);

    string dataset = Read::Read_String_from_File(age_dataset, "galaxy_age_dataset");
    data_info.GetInfo(dataset);

	string age_data = Read::Read_String_from_File(age_dataset, "galaxy_age_file");

	Read_Data( age_data );

	gft_marg	= Read::Read_Bool_from_File(age_dataset, "gft_marg");
}

//	Age data format: z   age   age_err

void Data_Age::Read_Data( std::string& data_file ){

	vector<double> vz;
	vector<double> vage;
	vector<double> vage_err;

	ifstream infile(data_file.c_str());

	if( !infile.good() ){
		std::string err = "\n*** Data_Age::Read_Data( std::string data_file ) ==> failed to open data file: " + data_file;
		throw runtime_error(err);
	}

	double zi, agei, age_erri;
	string line;

	age_num = 0;
	while( getline(infile, line) ){
		if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){

			stringstream ss(line);
			ss >> zi >> agei >> age_erri;

			vz.push_back(zi);
			vage.push_back(agei);
			vage_err.push_back(age_erri);

			age_num = age_num + 1;
		}
	}

	infile.close();

	if( age_num > 0 ){

		z		= new double[age_num];
		age		= new double[age_num];
		age_err	= new double[age_num];

		for(int i=0; i<age_num; ++i){
			z[i]		= vz[i];
			age[i]	= vage[i];
			age_err[i]= vage_err[i];
		}

		age_init = true;
	}
	else{
		string err = "\n*** Data_Age::Read_Data( std::string data_file ) ==> no age data is readed ... sorry,  stop here ...";
		throw runtime_error(err);
	}
}
