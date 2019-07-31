///////////////////////////////////////////////////////////////////////////////////////
//  This program is for testing ClassMC, for details pls see CosmologyTest.cpp
//
//  @ Sep-11-2017
//  Youhua Xu
///////////////////////////////////////////////////////////////////////////////////////

#include "ClassMC.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

CosmoTheory theory;
DataList    datalist;

// vector<string> mcmc_params;
// vector<string> derived_params;


int main( int argc, char* argv[] ){

    MPI::Init(argc, argv);

    int rank = MPI::COMM_WORLD.Get_rank();

    if( rank == 0 ){
        cout << "\n ====================\n";
        cout << " running ClassMC test";
        cout << "\n ====================\n";
    }

    if( argc < 2 ){
        if( rank == 0 ){
            cout << "\n"
                 << argv[0]
                 << ": no input parameter file found, exit now! \n";
        }

        MPI::Finalize();
		return 0;
    }

//  get configuration file
    string paramfile(argv[1]);

//  configure datasets.
	datalist.Init(paramfile);
	datalist.compute_sigma8 = true;

//  configure cosmological model
    theory.Init(paramfile);

//  decide what to calculate and pass cosmological MCMC parameters to datalist.ParList
    theory.Init(datalist);

//	Test part
    CosmoTheoryTest CTT;
    CTT.RunTest(paramfile,&theory);

    // theory.engine->computeCls();

//     int N = 1000;
//     double zmin = 1e-4;
//     double zmax = 3.0;
    
//     double lnzmin = log(zmin);
//     double lnzmax = log(zmax);

// //    ofstream outfile("classmc_z_Dl_mu_LCDM.txt");
// 	ofstream outfile("yyyyyyy.txt");
//     outfile << "# "
//             << setw(10) << "zi" << "\t"
//             << setw(15) << "Dl" << "\t"
//             << setw(15) << "mu" << endl;

//     for( int i=0; i<N; ++i ){
//         double lnzi = lnzmin + i*(lnzmax-lnzmin)/(N-1.);
//         double zi = exp(lnzi);
//         double Dl = theory.engine->get_Dl(zi);
//         double mu = 5.*log10(Dl) + 25.;
//         cout << setw(12) << setprecision(10) << zi << "\t"
//              << setw(15) << setprecision(10) << Dl << "\t"
//              << setw(15) << setprecision(10) << mu << endl;

//         setprecision(10);
//         outfile << setw(12) << setprecision(10) << zi << "\t"
//                 << setw(15) << setprecision(10) << Dl << "\t"
//                 << setw(15) << setprecision(10) << mu << endl;
//     }

//     outfile.close();

    MPI::Finalize();
}
