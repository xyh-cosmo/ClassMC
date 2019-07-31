#include <imcmc/parser++.hpp>
#include "CosmoTheory.hpp"
#include "Hubble.hpp"

using namespace std;
using namespace imcmc::parser;

Data_Hubble::Data_Hubble(){
	size        = 0;
	has_data    = false;
}

Data_Hubble::~Data_Hubble(){
	if( has_data ){
	    delete[] z;
	    delete[] Hz;
	    delete[] dHz;
	}
}

void Data_Hubble::ReadData( string& paramfile ){

	int my_rank = MPI::COMM_WORLD.Get_rank();

	string datafile;
	bool include_err = false;
	bool use_mock_hubble = false;

	if( Read::Has_Key_in_File(paramfile, "datafile") ){
	    datafile = Read::Read_String_from_File(paramfile, "datafile");
	}
    else{
        string err = "\n";
        err += "**** counl not find key: datafile !";
//        throw runtime_error(err);
        MPI::COMM_WORLD.Abort(MPI::COMM_WORLD.Get_rank());
    }

    data_info.GetInfo(datafile);

    if( Read::Has_Key_in_File(paramfile, "use_mock_hubble") ){

        use_mock_hubble = Read::Read_Bool_from_File(paramfile,"use_mock_hubble");

        if( use_mock_hubble == true ){
            if( Read::Has_Key_in_File(paramfile,"mock_hubble_has_error") )
                include_err = Read::Read_Bool_from_File(paramfile, "mock_hubble_has_error");
        }
    }

	vector<double> ztemp, Hztemp, dHztemp;
	ifstream infile( datafile.c_str() );

	string line;

	if( use_mock_hubble == true ){
	    while( getline(infile, line) ){
	        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
	            double zi, Hzi, dHzi, HzNoErri=0;
	            stringstream stream(line);
	            stream >> zi >> Hzi >> dHzi >> HzNoErri;
	            ztemp.push_back(zi);

				if( fabs(HzNoErri) < 1.0 ){
					string err = "\n";
					err += "==============================================\n";
					err += "> ****  error in mock hubble data file!  *****\n";
					err += "> you may need to set mock_hubble_has_error = T,\n";
					err += "> or do not use mock hubble data.\n";
					err += "==============================================\n";
					throw runtime_error(err);
				}

			    if( include_err )
	            	Hztemp.push_back(Hzi);
			    else
				    Hztemp.push_back(HzNoErri);

	            dHztemp.push_back(dHzi);
	            ++size;
	        }
	    }
	}
	else{
	    while( getline(infile, line) ){
	        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
	            double zi, Hzi, dHzi;
	            stringstream stream(line);
	            stream >> zi >> Hzi >> dHzi;
	            ztemp.push_back(zi);
            	Hztemp.push_back(Hzi);
	            dHztemp.push_back(dHzi);
	            ++size;
	        }
	    }
	}

	if( size > 1 ){
	    has_data = true;

	    z   = new double[size];
	    Hz  = new double[size];
	    dHz = new double[size];

	    for(int i=0; i<size; ++i){
	        z[i]    = ztemp[i];
	        Hz[i]   = Hztemp[i];
	        dHz[i]  = dHztemp[i];
	    }
	}

	//  NOTE: not all Hz data will be used, this depends on which redshfit range were chosen to test
	//  our theories or models, this was done by specifying a max-redshift value to key "hubble_z_max"

	if( Read::Has_Key_in_File(paramfile, "hubble_z_max") ){
	    if( Read::Has_Value(paramfile, "hubble_z_max", "double") ){
	        hubble_z_max = Read::Read_Double_from_File(paramfile, "hubble_z_max");
	    }
	    else{
	        MPI_cout("void Data_Hubble::ReadData() ==> found key \"hubble_z_max\" in the", my_rank);
	        MPI_cout("paramfile: " + paramfile, my_rank);
	        MPI_cout("but no max-z were given, so the hubble_z_max will be set to the max-z", my_rank);
	        MPI_cout("of the input data", my_rank);

	        hubble_z_max = z[0];
	        for(int i=1; i<size; ++i){
	            if( hubble_z_max <= z[i] )
	                hubble_z_max = z[i];
	        }
	    }
	}
	else{
	    MPI_cout("void Data_Hubble::ReadData() ==> you have not specify \"hubble_z_max\" in the", my_rank);
	    MPI_cout("paramfile: " + paramfile, my_rank);
	    MPI_cout("so the hubble_z_max will be set to the max-z of the input data", my_rank);

	    hubble_z_max = z[0];
	    for(int i=1; i<size; ++i){
	        if( hubble_z_max <= z[i] )
	            hubble_z_max = z[i];
	    }
	    cout << "==> Max z of H(z) data is " << hubble_z_max << endl;
	}

	infile.close();
}

void Data_Hubble::Init( string& paramfile ){
	ReadData( paramfile );
}
