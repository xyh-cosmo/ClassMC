#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "multinest.h"

#include "ClassMC.hpp"

#define NEG_MAX -1E90

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

CosmoTheory     theory;
DataList        datalist;

vector<string>  mcmc_params;
vector<string>  derived_params;

vector<string>  output_mcmc_params;
vector<string>  output_derived_params;

imcmc_double    all_params;

vector<string>  mcmc_params_float;  //  these are the true parameters to be sampled
imcmc_double    mcmc_params_lower;
imcmc_double    mcmc_params_upper;

istate          state;


/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
    string par;
    double chisq,lndet;
    double x_lower=0, x_upper=0;

    for( int i=0; i<ndim; ++i ){    // update mcmc_parameter from *Cube
        par = mcmc_params_float[i];
        x_lower = mcmc_params_lower[par];
        x_upper = mcmc_params_upper[par];
        all_params[par] = x_lower + (x_upper-x_lower)*Cube[i];
//        Cube[i] = all_params[par];  // update MCMC parameters for output
    }

    // call the total likelihood function & update derived parameters
    Likelihoods(all_params,
                lndet,
                chisq,
                &theory,
                &datalist,
                state);

    if( state.this_like_is_ok == false ){
        lnew = NEG_MAX*10;
        return;
    }

    lnew = -0.5*chisq;

    if( output_mcmc_params.size() == npars ){
        for( int i=0; i<ndim; ++i ){
            Cube[i] = all_params[output_mcmc_params[i]];  // update MCMC parameters for output
        }
    } else {
        cout << "Error in LogLike: output_mcmc_params.size() != ndim\n";
        MPI::Finalize();
        exit(0);
    }

    // update derived parameters
//    if( derived_params.size() > 0 ){
//        for( int i=ndim; i<npars; ++i ){
//            par = derived_params[i-ndim];
//            Cube[i] = all_params[par];
//        }
//    }
}

/***********************************************************************************************************************/




/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays


	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns

	int i, j;

	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];



	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column

	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{

    MPI::Init(argc, argv);

    int rank = MPI::COMM_WORLD.Get_rank();

    if( argc < 2 ){

        if( rank == 0 ){
            cout << "====================================\n";
            cout << "=> No Input Parameter File Found !\n";
            cout << "=> So stop here !\n";
            cout << "====================================\n";
        }

        MPI::Finalize();
        return 0;
    }

//      ============================================================================

//  get configuration file:
    string paramfile(argv[1]);

//  configure datasets.
    datalist.Init(paramfile);

//  configure cosmological model
    theory.Init(paramfile);

//  decide what to calculate and pass cosmological MCMC parameters to datalist.ParList
    theory.Init(datalist);

    datalist.Get_MCMC_Params(mcmc_params);
    datalist.Get_Derived_Params(derived_params);

//  read lower and upper bounds
    for( int i=0; i<mcmc_params.size(); ++i ){
        int nvals=-1;
        double *vals = NULL;
        vals = Read::Read_Array_of_Double_from_File(paramfile, mcmc_params[i], nvals );
//        cout << "@debug --> nvals = " << nvals << endl;
        if( nvals == 3 ){
            mcmc_params_float.push_back(mcmc_params[i]);
            mcmc_params_lower[mcmc_params[i]] = vals[1];
            mcmc_params_upper[mcmc_params[i]] = vals[2];

            if( rank == 0 ){
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


    int mcmc_par_num = mcmc_params_float.size();
    int derived_par_num = derived_params.size();

    //  configure orders of output params:
    if( Read::Has_Key_in_File( paramfile, "output_params" ) ) {
        int nvalue = Read::Num_of_Value_for_Key( paramfile, "output_params" );
        std::string *name = new std::string[nvalue];
        Read::Read_Array_from_File(paramfile, "output_params", name, nvalue);
        for( int i=0; i<nvalue; ++i ){
           output_mcmc_params.push_back(name[i]);
        }
    } else {
        cout << "cannot found keyword \'output_params\' in : " << paramfile << endl;
        MPI::Finalize();
        exit(0);
    }

    //  configure orders of output params:
//    if( Read::Has_Key_in_File( paramfile, "output_dparams" ) ) {
//        int nvalue = Read::Num_of_Value_for_Key( paramfile, "output_dparams" );
//        std::string *name = new std::string[nvalue];
//        Read::Read_Array_from_File(paramfile, "output_dparams", name, nvalue);
//        for( int i=0; i<nvalue; ++i ){
//           output_derived_params.push_back(name[i]);
//        }
//    }

    for( int i=0; i<mcmc_params_float.size(); ++i ){
        if( rank == 0 )
            cout << "debug: found MCMC par = " << mcmc_params_float[i] << endl;
    }

    for( int i=0; i<derived_params.size(); ++i ){
        if( rank == 0 ){
            cout << "debug: found derived par = " << derived_params[i] << endl;
        }

        all_params[derived_params[i]] = -999;   // add derived params into all_params, -999 is just a place holder ...
        output_mcmc_params.push_back(derived_params[i]);
    }

//    exit(0);

//      ============================================================================

	// set the MultiNest sampling parameters


	int IS = 1;					// do Nested Importance Sampling?

	int mmodal = 0;					// do mode separation?

	int ceff = 1;					// run in constant efficiency mode?

	int nlive = 1000;				// number of live points

	double efr = 0.8;				// set the required efficiency

	double tol = 0.001;				// tol, defines the stopping criteria

	int ndims = mcmc_par_num;					// dimensionality (no. of free parameters)

	int nPar = mcmc_par_num + derived_par_num;					// total no. of parameters including free & derived parameters

	int nClsPar = 0;				// no. of parameters to do mode separation on

	int updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations

	double Ztol = NEG_MAX;				// all the modes with logZ < Ztol are ignored

	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)

	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;

	char root[1000]; // = "chains/eos_cdm";			// root for output files

	string root_str = Read::Read_String_from_File(paramfile,"chain_root");
	sprintf(root,"%s",root_str.c_str());

	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock

	int fb = 1;					// need feedback on standard output?

	int resume = 0;					// resume from a previous job?

	int outfile = 1;				// write output files?

	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization

	double logZero = NEG_MAX;				// points with loglike < logZero will be ignored by MultiNest

	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

	void *context = 0;				// not required by MultiNest, any additional information user wants to pass



	// calling MultiNest

	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLike, dumper, context);

	MPI::Finalize();
}

/***********************************************************************************************************************/
