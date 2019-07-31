#include "RSD.hpp"
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>

using namespace std;
using namespace imcmc::parser;

Data_RSD::Data_RSD(){
    int size = -1;
    zeff = NULL;
    val  = NULL;
    err  = NULL;
}

Data_RSD::~Data_RSD(){
    if( zeff !=NULL ){
        delete[] zeff;
        zeff = NULL;
    }

    if( val != NULL ){
        delete[] val;
        val = NULL;
    }

    if( err != NULL ){
        delete[] err;
        err = NULL;
    }
}

int Data_RSD::Init( std::string& rsd_dataset ){
    data_info.GetInfo(rsd_dataset);
    string datafile = Read::Read_String_from_File(rsd_dataset,"data_file");
    Read_Data(datafile);
    return 0;
}

int Data_RSD::Read_Data( std::string& datafile ){
    int data_size = 0;

    vector<double> zeff_tmp, val_tmp, err_tmp;
    string line;
    ifstream infile(datafile.c_str());
    while( getline(infile,line) ){
        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
            double zi, vali, erri;
            stringstream stream(line);
            stream >> zi >> vali >> erri;
            zeff_tmp.push_back(zi);
            val_tmp.push_back(vali);
            err_tmp.push_back(erri);
            ++data_size;

            // cout << "z = " << zi << ", val = " << vali << ", err = " << erri << endl;
        }
    }

    // cout << "total # of rsd data is: " << data_size << endl;
    // exit(0);

    infile.close();

    size = data_size;
    zeff = new double[size];
    val  = new double[size];
    err  = new double[size];

    for( int i=0; i<data_size; ++i ){
        zeff[i] = zeff_tmp[i];
        val[i] = val_tmp[i];
        err[i] = err_tmp[i];
    }

    return data_size;
}

