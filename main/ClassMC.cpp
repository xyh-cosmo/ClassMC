//  ==================================================================================
//  Some notes:
//  CosmoTheory has to be initialized in the following order:
//  1) theory.init(argv[1])
//  2) theory.init(&data_plk)
//  3) theory.init()
//
//  step 1) is to copy all relevant cosmological parameters into CosmoParams pars.
//  step 2) is to determine what Cls are going to be computed, according user provided
//          Planck data combinations.
//  step 3) is the real step of creating a CosmoTheory object, which is actually ClassEngine
//  ==================================================================================

#include "ClassMC.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

CosmoTheory theory;
DataList    datalist;

vector<string> mcmc_params;
vector<string> derived_params;


int main( int argc, char* argv[] ){

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

    ensemble_workspace ew;

//  get configuration file:
    string paramfile(argv[1]);

//  configure datasets.
	datalist.Init(paramfile);

//  PLANCK data has a self-check, which needs some time. We need to make sure
//  all datasets are properly initialized before we proceed
    MPI::COMM_WORLD.Barrier();

//  configure cosmological model
    theory.Init(paramfile);

//  decide what to calculate and pass cosmological MCMC parameters to datalist.ParList
    theory.Init(datalist);

    datalist.Get_MCMC_Params(mcmc_params);

    datalist.Get_Derived_Params(derived_params);

	MPI::COMM_WORLD.Barrier();

    ew.add_likelihood(  Likelihoods,
                        mcmc_params,
                        derived_params,
                        &theory,
                        &datalist );

    ew.init(paramfile);

    ew.do_sampling();

    MPI::Finalize();
}
