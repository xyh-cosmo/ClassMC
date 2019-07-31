#include "BAO.hpp"
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>

using namespace std;
using namespace imcmc::parser;

Data_BAO_Base::Data_BAO_Base(){
    type = -1;
    size = -1;
    z   = NULL;
    val = NULL;
    err = NULL;
    icov= NULL;
}

Data_BAO_Base::~Data_BAO_Base(){

    if( z != NULL ){
        delete[] z;
        z = NULL;
    }

    if( val != NULL ){
        delete[] val;
        val = NULL;
    }

    if( err != NULL ){
        delete[] err;
        err = NULL;
    }

    if( (size > 1) && (icov != NULL) ){
        delete[] icov;
        icov = NULL;
    }

}


int Data_BAO_Base::read_data( string& datafile ){

    // cout << "==> read BAO data\n";

    ifstream data(datafile.c_str());
    string line;

    int bao_size = 0;
    vector<double> bao_z;
    vector<double> bao_val;
    vector<double> bao_err;

    while( getline(data,line) ){
        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
            stringstream stream(line);
            double zi,vi,ei;
            stream >> zi >> vi >> ei;
            bao_z.push_back(zi);
            bao_val.push_back(vi);
            bao_err.push_back(ei);
            ++bao_size;
        }
    }

    data.close();

    if( bao_size > 0 ){
        size = bao_size;
        z   = new double[size];
        val = new double[size];
        err = new double[size];
        // icov is to be handeled in another member function
        for( int i=0; i<size; ++i ){
            z[i]    = bao_z[i];
            val[i]  = bao_val[i];
            err[i]  = bao_err[i];
        }
    }
    else{
		string errmsg = "found no BAO data in datafile: " + datafile;
        MPI_ClassMC_ERROR(errmsg);
    }

    return bao_size;
}

// if this member function is called, then we must have size>0
int Data_BAO_Base::read_icov( string& icovfile ){

    // cout << "==> read BAO cov\n";

    this->icov = new double[size*size];

    ifstream fp_icov(icovfile.c_str());
    string line;
    
    if( fp_icov.good() == false ){
	    string errmsg = "found no BAO (inv-)covariance data in icovfile: " + icovfile;
        MPI_ClassMC_ERROR(errmsg);
    }

    int cnt_line=0;
    while( getline(fp_icov,line) ){
        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){

            int cnt=0;
            double temp;
            vector<double> icovij;
            stringstream stream(line);

            while( stream >> temp ){
                icovij.push_back(temp);
                ++cnt;
            }

            if( cnt == size ){
                for( int j=0; j<size; ++j ){
                    this->icov[cnt_line*size+j] = icovij[j];
                }
            }
            else{
                cout << "\n==> wrong number of icov values in line: " << cnt_line+1 << endl;
                MPI::COMM_WORLD.Abort(MPI::COMM_WORLD.Get_rank());
            }

            ++cnt_line;
        }        
    }

    fp_icov.close();

    return size;
}


Data_BAO::Data_BAO(){
    bao_data_size = 0;
    use_Hu_fitting = true;
}

Data_BAO::~Data_BAO(){

    for( size_t i=0; i != BAO_List.size(); ++i ){
        BAO_List[i]->~Data_BAO_Base();
        BAO_List[i] = NULL;
    }
}

void Data_BAO::Init( std::string& bao_dataset ){

    data_info.GetInfo(bao_dataset);

    int max_bao_index=-1;
    if( Read::Has_Key_in_File(bao_dataset,"max_bao_index") ){
        max_bao_index = Read::Read_Int_from_File(bao_dataset,"max_bao_index");
        if( max_bao_index <= 0 ){
            cout << "\n==> max_bao_index should not be less than zero.\n";
            MPI::COMM_WORLD.Abort(MPI::COMM_WORLD.Get_rank());
        }
    }

    for( int i=1; i<=max_bao_index; ++i ){

        string idx = Read::IntToString(i);

        string fname1 = "data_file_"+idx;

        int data_size = -1; // -1 can be used to check data
        if( Read::Has_Key_in_File(bao_dataset,fname1) ){

            Data_BAO_Base *bao_base;
            bao_base = new Data_BAO_Base;

            string datafile = Read::Read_String_from_File(bao_dataset,fname1);
            data_size = bao_base->read_data(datafile);

            if( data_size <= 0 ){
                cout << "\n==> " + datafile << " contains no BAO data, so skip it!\n";
                bao_base->~Data_BAO_Base();
                bao_base = NULL;
            }
            else{

                // ok, we got some data, now get its type
                if( Read::Has_Key_in_File(bao_dataset,"data_type_"+idx)){
                    bao_base->type = Read::Read_Int_from_File(bao_dataset,"data_type_"+idx);
                }
                else{
                    cout << "\n==> failed to get bao data_type for: " << fname1 << endl;
                    MPI::COMM_WORLD.Abort(MPI::COMM_WORLD.Get_rank());
                }

                if( data_size > 1 ){
                    string fname2 = "data_icov_"+idx;
                    if( Read::Has_Key_in_File(bao_dataset,fname2) ){
                        string icovfile = Read::Read_String_from_File(bao_dataset,fname2);
                        int icov_size = bao_base->read_icov(icovfile);
                        if( icov_size != data_size ){
                            bao_base->~Data_BAO_Base();
                            cout << "\n==> readed BAO inverse covariance matrix has wrong size !\n";
                            MPI::COMM_WORLD.Abort(MPI::COMM_WORLD.Get_rank());
                        }
                    }
                }

                // trying to get BAO data name
                string data_name;
                if( Read::Has_Key_in_File(bao_dataset,"data_name_"+idx)){
                    data_name = Read::Read_String_from_File(bao_dataset,"data_name_"+idx);
                    data_info.DataName += (" ,"+data_name);
                }
                else{
                    cout << "\n==> failed to get bao data data_name for: " << fname1 << endl;
                    MPI::COMM_WORLD.Abort(MPI::COMM_WORLD.Get_rank());
                }

                BAO_List.push_back(bao_base);
                ++bao_data_size;
            }
        }

    }

    if( Read::Has_Key_in_File(bao_dataset,"use_Hu_fitting") ){
        use_Hu_fitting = Read::Read_Bool_from_File(bao_dataset,"use_Hu_fitting");
    }
    else{
        MPI_ClassMC_WARNING("use_Hu_fitting is not found in:"+bao_dataset+", so set it to default value \'true\'");
    }

}


