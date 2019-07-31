//	========================================================================
//	Author: Youhua Xu
//	Date:	Jan-7, 2017
//
// 	Structure defined in this file is used to handel ALL parameters that
// 	related to Cosmology, including both MCMC parameters and derived ones.
//	========================================================================

#ifndef __PARAMLIST__
#define __PARAMLIST__

#include <vector>
#include <string>
#include <set>
#include <map>

/*
    Store information about a specific data set
*/
struct DataInfo{ // use name DataParam will be better ...

    int rank;

    std::string DataName;
    std::map<std::string,std::string> MCMC_Params;
    std::map<std::string,std::string> Derived_Params;

    DataInfo();
    ~DataInfo();
    void GetInfo( std::string& dataset );

//  if some dataset needs some specific parameters, for example, CLik needs
//  different extra parameters when different likelihood files are used.  In this
//  case we do not specify the MCMC & derived parameters inside the dataset.ini,
//  but directly add those parameters into DataInfo using the following two methods:
    void Add_MCMC_Param( std::string& mcmc_param );
    void Add_Derived_Param( std::string& derived_params );
};

//  sometimes, one dataset may not be used together with another one, i.e., RSD vs. BAO ...
struct DataConflict{
//  NOT implemented yet ..
};

/*
    Collect ALL mcmc & derived parameters
*/
struct ParamList {

    int rank;

    std::string my_name;

    std::map<std::string, std::string> MCMC_Params;
    std::map<std::string, std::string> Derived_Params;

    ParamList(std::string& paramlistname );
    ParamList();
    ~ParamList();

//	Member functions: the first two are for
    int Add_MCMC_Params( std::map<std::string,std::string>& mcmc_params,std::string& datasetname );
    int Add_Derived_Params( std::map<std::string,std::string>& derived_params,std::string& datasetname );

//  add params from DataInfo strcuture(s)
    int Add_Params( DataInfo& datainfo );

//  these two are to used inside NON-data structures
    int Add_MCMC_Param(std::string& mcmc_param);    //  add a single MCMC param
    int Add_Derived_Param(std::string& derived_param);

//	For debug
    void Print_MCMC_Params();
    void Print_Derived_Params();

    void Get_MCMC_Params( std::vector<std::string>& mcmc_params );
    void Get_Derived_Params( std::vector<std::string>& derived_params );
};

#endif //__PARAMLIST__
