#include <stdio.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "multinest.h"
#include "CC_ini_likelihood.hpp"
#include "interfaces.hpp"
#include "ClassMC.hpp"

#define NEG_MAX -1E90

#define MPI_exit {MPI::Finalize(); exit(0);}

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

/*****************************************************************/
/*                          Global vaiables                      */
CosmoTheory     theory;
DataList        datalist;

imcmc_double    all_params;
imcmc_double    mcmc_params_lower;  //  store upper bounds of MCMC parameters
imcmc_double    mcmc_params_upper;  //  store lower bounds of MCMC parameters

vector<string>  mcmc_params;        //  store all MCMC parameters needed by the model and likelihood, but some of them might be fixed to a given value
vector<string>  derived_params;     //  store derived parameters

vector<string>  mcmc_params_float;  //  these are the true parameters to be sampled
vector<string>  output_params;      // contains only the float MCMC parameters; the derived parameters are not handled here!!!!
vector<string>  sorted_mcmc_params; // store MCMC parameter names acoording to the order specified by USER !

Settings settings;

/*****************************************************************/
double loglikelihood (double theta[], int nDims, double phi[], int nDerived)
{
    double logL=0;
    double lndet=0,chisq=0;

    for( int i=0; i<nDims; ++i ){
        all_params[sorted_mcmc_params[i]] = theta[i];
    }

    istate state;    
    Likelihoods(all_params,lndet,chisq,&theory,&datalist,state);

    if( state.this_like_is_ok == false ){
        logL = NEG_MAX*10;
        return logL;
    }

    // update derived parameters (these parameters are actually updated inside "theory")
    for( int i=0; i<nDerived; ++i ){
        phi[i] = all_params[derived_params[i]];
    }
    
    logL = -0.5*chisq;
    
    return logL;

}

// Prior function
void prior (double cube[], double theta[], int nDims)
{
    for(int i=0;i<nDims;i++){
        string par = sorted_mcmc_params[i];
        double x_lower = mcmc_params_lower[par];
        double x_upper = mcmc_params_upper[par];
        theta[i] = x_lower + (x_upper-x_lower)*cube[i];
    }
}

// Dumper function
void dumper(int ndead,int nlive,int npars,double* live,double* dead,double* logweights,double logZ, double logZerr)
{
//  do nothing
}


// Ini path reading function
//
// If you want the file path to the ini file (for example for pointing to other config files), store ite value in a global variable with this function
void set_ini(std::string ini_str_in)
{
    // do nothing
}

// Setup of the loglikelihood
void setup_loglikelihood()
{
}

//============================================================
int main(int argc, char *argv[])
{
    MPI::Init(argc, argv);

    string paramfile = argv[1];     //  configuration filename
    datalist.Init(paramfile);       //  configure datasets.
    theory.Init(paramfile);         //  configure cosmological model
    theory.Init(datalist);          //  decide what to calculate and pass cosmological MCMC parameters to datalist.ParList
    datalist.Get_MCMC_Params(mcmc_params);
    datalist.Get_Derived_Params(derived_params);
    
//  read lower and upper bounds of the MCMC parameters
    for( int i=0; i<mcmc_params.size(); ++i ){
        int nvals=-1;
        double *vals = NULL;
        vals = Read::Read_Array_of_Double_from_File(paramfile, mcmc_params[i], nvals ); 

        if( nvals == 3 ){
            mcmc_params_float.push_back(mcmc_params[i]);
            mcmc_params_lower[mcmc_params[i]] = vals[1];
            mcmc_params_upper[mcmc_params[i]] = vals[2];
            
            if( MPI::COMM_WORLD.Get_rank() == 0 ){
                cout << " Par: " << setw(20) << mcmc_params[i]
                     << "\tlower = " << setw(10) << mcmc_params_lower[mcmc_params[i]] 
                     << "\tupper = " << setw(10) << mcmc_params_upper[mcmc_params[i]] << endl;
            }
        }
        
        all_params[mcmc_params[i]] = vals[0];
        
        if( vals != NULL ){
            delete[] vals;
        }
    }

    //  configure orders of output params:
    if( Read::Has_Key_in_File( paramfile, "output_params" ) ) {
        int nvalue = Read::Num_of_Value_for_Key( paramfile, "output_params" );
        
        if( nvalue <= 0 ){
            if( MPI::COMM_WORLD.Get_rank() == 0 ){
                cout << "Error: parameters to be output should be greater than zero!\n";
            }
            MPI_exit;
        }
        
        std::string *name = new std::string[nvalue];
        Read::Read_Array_from_File(paramfile, "output_params", name, nvalue);
        for( int i=0; i<nvalue; ++i ){
           output_params.push_back(name[i]);
        }
    } else {
        if( MPI::COMM_WORLD.Get_rank() == 0 ){
            cout << "cannot found keyword \'output_params\' in : " << paramfile << endl;
        }
        MPI_exit;
    }
    
    if( output_params.size() > mcmc_params_float.size() ){  //  ready to sort the MCMC parameters:
        if( MPI::COMM_WORLD.Get_rank() == 0 ){
            cout << "Error: output_params.size() should NOT be greater than mcmc_params_float.size(), stop!\n";
        }
        MPI_exit;
    }else if( output_params.size() == mcmc_params_float.size() ){
        for( int i=0; i<output_params.size(); ++i ){
        //  check if output_params[i] is in mcmc_params_float[i];
            vector<string>::iterator it;
            it = find(mcmc_params_float.begin(),mcmc_params_float.end(),output_params[i]);
            if( it != mcmc_params_float.end() ){
                sorted_mcmc_params.push_back(output_params[i]);
            } else {
                if( MPI::COMM_WORLD.Get_rank() == 0 ){
                    cout << "Error: " << output_params[i] << " is not in mcmc_params_float!\n";
                }
                MPI_exit;
            }
        }
    }else if( output_params.size() < mcmc_params_float.size() ){
        for( int i=0; i<output_params.size(); ++i ){
        //  check if output_params[i] is in mcmc_params_float
            vector<string>::iterator it;
            it = find( mcmc_params_float.begin(), mcmc_params_float.end(), output_params[i] );
            if( it != mcmc_params_float.end() ){
                sorted_mcmc_params.push_back (output_params[i] );
            } else {
                if( MPI::COMM_WORLD.Get_rank() == 0 ){
                    cout << "Error: " << output_params[i] << " is not in mcmc_params_float!\n";
                }
                MPI_exit;
            }
        }
        
        for( int i=0; i<mcmc_params_float.size(); ++i ){    //  add the rest MCMC parameters ...
        //  check if mcmc_params_float[i] is in output_params
            vector<string>::iterator it;
            it = find( output_params.begin(), output_params.end(), mcmc_params_float[i]);
            if( it == output_params.end() ){
                sorted_mcmc_params.push_back(mcmc_params_float[i]);
            }
        }
        
        //  check if sorted_mcmc_params.size() == mcmc_params_float.size()
        if( sorted_mcmc_params.size() != mcmc_params_float.size() ){
            if( MPI::COMM_WORLD.Get_rank() == 0 ){
                cout << "Fatal Error: sorted_mcmc_params.size() != mcmc_params_float.size()\n";
            }
            MPI_exit;
        }
    }
    
    for( int i=0; i<sorted_mcmc_params.size(); ++i ){
        if( MPI::COMM_WORLD.Get_rank() == 0 ){
            cout << "==> sorted MCMC parameter : " << sorted_mcmc_params[i] << endl;
        }
    }

    for( int i=0; i<derived_params.size(); ++i ){
        if( MPI::COMM_WORLD.Get_rank() == 0 ){
            cout << "==> derived parameter : " << derived_params[i] << endl;
        }
        
        all_params[derived_params[i]] = -999;   // add derived params into all_params, -999 is just a place holder ...
    }
    
//    MPI_exit;

    int nDims, nDerived;
    nDims = mcmc_params_float.size();
    nDerived = derived_params.size();

    Settings settings(nDims,nDerived);

    settings.nlive         = 500;
    settings.num_repeats   = settings.nDims*2;
    settings.do_clustering = false;

    settings.precision_criterion = 1e-3;
    settings.logzero        = -1e30;

    settings.base_dir      = "chains";
    settings.file_root     = "test";

    settings.write_resume  = false;
    settings.read_resume   = false;
    settings.write_live    = true;
    settings.write_dead    = false;
    settings.write_stats   = false;

    settings.equals        = false;
    settings.posteriors    = true;
    settings.cluster_posteriors = false;

    settings.feedback      = 1;

//    settings.compression_factor  = 0.36787944117144233;
//    settings.boost_posterior= 10.0;

//    setup_loglikelihood();
    run_polychord(loglikelihood,prior,dumper,settings) ;

    MPI::Finalize();
}

