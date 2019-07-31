#include <iostream>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "mpi.h"

#include "ParamList.hpp"
using namespace std;
using namespace imcmc::parser;

//  ===============
//  struct DataInfo
//  ===============

DataInfo::DataInfo(){
    rank = MPI::COMM_WORLD.Get_rank();
    DataName = "Unknown";
    MCMC_Params.clear();
    Derived_Params.clear();
}

DataInfo::~DataInfo(){
    MCMC_Params.clear();
    Derived_Params.clear();
}

void DataInfo::GetInfo(string& dataset ){

    if( rank == ROOT_RANK )
        cout << "==> Getting data information from: " << dataset << endl;

    if( Read::Has_Key_in_File(dataset, "DataName") ){
        DataName = Read::Read_String_from_File(dataset,"DataName");
    }

    if( Read::Has_Key_in_File(dataset,"MCMC_Params") ){
        int nval = Read::Num_of_Value_for_Key(dataset, "MCMC_Params");
        if( nval > 0 ){
            string val;
            string *name = new string[nval];
            Read::Read_Array_of_String_from_File( dataset, "MCMC_Params", name, nval );
            for( int i=0; i<nval; ++i ){
                Read::Read_Value_from_File(dataset,name[i],val);
                MCMC_Params[name[i]] = val;
            }
        }
    }

    if( Read::Has_Key_in_File(dataset,"Derived_Params") ){
        int nval = Read::Num_of_Value_for_Key(dataset, "Derived_Params");
        if( nval > 0 ){
            string val;
            string *name = new string[nval];
            Read::Read_Array_of_String_from_File( dataset, "Derived_Params", name, nval );
            for( int i=0; i<nval; ++i ){
                Read::Read_Value_from_File(dataset,name[i],val);
                Derived_Params[name[i]] = val;
            }
        }
    }
}

void DataInfo::Add_MCMC_Param( string& mcmc_param ){
    if( MCMC_Params.count(mcmc_param) == 0 ){
        MCMC_Params[mcmc_param] = "Unknown value";
    }
}

void DataInfo::Add_Derived_Param( string& derived_param ){
    if( Derived_Params.count(derived_param) == 0 ){
        Derived_Params[derived_param] = "Unknown value";
    }
}

//  ================
//  struct ParamList
//  ================

ParamList::ParamList(string& paramlistname ) {
    rank = MPI::COMM_WORLD.Get_rank();
    my_name = paramlistname;
    MCMC_Params.clear();
    Derived_Params.clear();
}

ParamList::ParamList() {
    my_name = "ParamList_for_ClassMC";
    MCMC_Params.clear();
    Derived_Params.clear();
}

ParamList::~ParamList() {
    MCMC_Params.clear();
    Derived_Params.clear();
}

int ParamList::Add_MCMC_Params( map<string,string>& mcmc_params, string& datasetname ) {

    int cnt=0;

    if( mcmc_params.size() > 0 ){

        map<string,string>::iterator it;

        for( it = mcmc_params.begin(); it != mcmc_params.end(); ++it ) {
            if( MCMC_Params.size() == 0 || MCMC_Params.count(it->first) == 0 ) {
                MCMC_Params[it->first] = it->second;
                ++cnt;
            }
            else {
                if( rank == ROOT_RANK )
                    cout << "\n==> MCMC Param: " << it->first << " already exist.\n";
            }
        }

        if( rank == ROOT_RANK )
            cout << "==> added " << cnt << " MCMC parameters from Data: " << datasetname << "\n";
    }

    return cnt;
}

int ParamList::Add_Derived_Params( map<string,string>& derived_params, string& datasetname ) {

    int cnt=0;

    if( derived_params.size() > 0 ){

        map<string,string>::iterator it;
        for( it = derived_params.begin(); it != derived_params.end(); ++it ) {
            if( Derived_Params.count(it->first) == 0 ) {
                Derived_Params[it->first] = it->second;
                ++cnt;
            }
            else {
                if( rank == ROOT_RANK )
                    cout << "\n==> Derived_Param: " << it->first << " already exist.\n";
            }
        }

        if( rank == ROOT_RANK )
            cout << "==> added " << cnt << " Derived parameters from Data: " << datasetname << "\n";
    }

    return cnt;
}

int ParamList::Add_Params( DataInfo& datainfo ){

    int n1 = Add_MCMC_Params(datainfo.MCMC_Params, datainfo.DataName);
    int n2 = Add_Derived_Params(datainfo.Derived_Params, datainfo.DataName);

    return ( n1 + n2 );
}


int ParamList::Add_MCMC_Param( string& mcmc_param ){
    if( MCMC_Params.count(mcmc_param) == 0 )
        MCMC_Params[mcmc_param] = "Unknown value";

    return 0;
}


int ParamList::Add_Derived_Param( string& derived_param ){
    if( Derived_Params.count(derived_param) == 0 )
        Derived_Params[derived_param] = "Unknown value";

    return 0;
}


void ParamList::Print_MCMC_Params() {
    if( MCMC_Params.size() > 0 ) {

        if( rank == ROOT_RANK )
            cout << "\n==> MCMC parameters stored in ParamList: " << my_name << " are:\n";

        map<string,string>::iterator it;
        for( it = MCMC_Params.begin(); it != MCMC_Params.end(); ++it ) {
            cout << setw(25) << it->first << " = " << it->second << "\n";
        }
    }
}

void ParamList::Print_Derived_Params() {
    if( Derived_Params.size() > 0 ) {

        if( rank == ROOT_RANK )
            cout << "\n==> Derived parameters stored in ParamList: " << my_name << " are:\n";

        map<string,string>::iterator it;
        for( it = Derived_Params.begin(); it != Derived_Params.end(); ++it ) {
            cout << setw(25) << it->first << " = " << it->second << "\n";
        }
    }
}

void ParamList::Get_MCMC_Params( std::vector<std::string>& mcmc_params ){
    if( mcmc_params.size() > 0 )
        mcmc_params.clear();

    map<string,string>::iterator it;
    for( it = MCMC_Params.begin(); it != MCMC_Params.end(); ++it ){
        mcmc_params.push_back(it->first);
    }
}

void ParamList::Get_Derived_Params( std::vector<std::string>& derived_params ){
    if( derived_params.size() > 0 )
        derived_params.clear();

    map<string,string>::iterator it;;
    for( it = Derived_Params.begin(); it != Derived_Params.end(); ++it ){
        derived_params.push_back(it->first);
    }
}
