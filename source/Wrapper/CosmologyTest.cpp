#include "ClassMC.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
// using namespace arma;

void CosmoTheoryTest::Init(string& paramfile){

    test_mu = false;
    test_BAO= false;
    test_Hz = false;
    test_weff= false;

    if( Read::Has_Key_in_File(paramfile,"test_mu") ){
        test_mu = Read::Read_Bool_from_File(paramfile,"test_mu");
    }

    if( Read::Has_Key_in_File(paramfile,"test_BAO") ){
        test_BAO = Read::Read_Bool_from_File(paramfile,"test_BAO");
    }

    if( Read::Has_Key_in_File(paramfile,"test_Hz") ){
        test_Hz = Read::Read_Bool_from_File(paramfile,"test_Hz");
    }

    if( Read::Has_Key_in_File(paramfile,"test_weff") ){
        test_weff = Read::Read_Bool_from_File(paramfile,"test_weff");
    }

    zmin = 1e-3;
    if( Read::Has_Key_in_File(paramfile,"test_zmin") ){
        zmin = Read::Read_Double_from_File(paramfile,"test_zmin");
    }

    zmax = 3.0; // this should be enough
    if( Read::Has_Key_in_File(paramfile,"test_zmax") ){
        zmax = Read::Read_Double_from_File(paramfile,"test_zmax");
    }

    output_size = 50;
    if( Read::Has_Key_in_File(paramfile,"test_output_size") ){
        output_size = Read::Read_Int_from_File(paramfile,"test_output_size");
    }

    test_root = "classmc_test";
    if( Read::Has_Key_in_File(paramfile,"test_output_root") ){
        test_root = Read::Read_String_from_File(paramfile,"test_output_root");
    }

    // test_root = "test_output/"+test_root;
}

void CosmoTheoryTest::RunTest(string& paramfile, CosmoTheory *cosmo ){

    Init(paramfile);

//	initialize redshifts where testing-quantities are to be evaluated
    double *z = new double[output_size];
    for( int i=0; i<output_size; ++i ){
        z[i] = zmin + (zmax-zmin)/(output_size-1.0)*i;
        double tmp = (log(1.+zmax)-log(1.+zmin))/(output_size-1.0)*i + log(1.+zmin);
        z[i] = exp(tmp)-1;
    }

    int cnt = 0;
    string fname;

    ofstream outfile;

//	testing distance modulus computation
    if( test_mu ){

        fname = test_root+"_mu.txt";
        outfile.open(fname.c_str());

        if( Read::Has_Key_in_File(paramfile,"SN_z_file") == true ){
            int z_num=0;
            vector<double> z;
            string SN_z_file = Read::Read_String_from_File(paramfile,"SN_z_file");
            ifstream infile(SN_z_file.c_str());
            string line;
            while( getline(infile,line) ){
                if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
                    double ztmp;
                    stringstream stream(line);
                    stream >> ztmp;
                    z.push_back(ztmp);
                    z_num++;
                }
            }
            
            int i=0;
            while( i < z_num ){
                double Dl = cosmo->engine->get_Dl(z[i]);
                outfile << setw(15) << setprecision(10) << z[i] << " "
                        << setw(15) << setprecision(10) << Dl << " "
                        << setw(15) << setprecision(10) << 5*log10(Dl)+25 << "\n";
                ++i;
            }
        }
        else{
            cout << "Cannot find key: mu_z_file in " << paramfile << endl;
            int i=0;
            while( i < output_size ){
                double Dl = cosmo->engine->get_Dl(z[i]);
                outfile << setw(15) << setprecision(10) << z[i] << " "
                        << setw(15) << setprecision(10) << Dl << " "
                        << setw(15) << setprecision(10) << 5*log10(Dl)+25 << "\n";
                ++i;
            }
        }
        outfile.close();
    }

    // exit(0);

//	testing BAO observables' calculation
    cnt = 0;
    if( test_BAO ){
        double zd,zr,rs,Dv,Om,Da;
        fname = test_root+"_BAO.txt";
        outfile.open(fname.c_str());
        outfile << "# z zrec zdrag rs  Om  Dv  Da\n";
        while( cnt < output_size ){
            zr = cosmo->engine->z_rec();

            // zd = cosmo->engine->z_drag();
            zd = cosmo->engine->z_drag_Hu();
            
            // rs = cosmo->engine->rs_drag();
            rs = cosmo->engine->get_rs(zd);
            
            Om = cosmo->engine->get_Omega_m();
            Dv = cosmo->engine->get_Dv(z[cnt]);
            Da = cosmo->engine->get_Da(z[cnt]);
            outfile << z[cnt] << "\t"
                    << zr << "\t"
                    << zd << "\t"
                    << rs << "\t"
                    << Om << "\t"
                    << Dv << "\t"
                    << Da << "\n";
            ++cnt;
        }
        outfile.close();
    }

//	testing Hublle parameter calculation
    cnt = 0;
    if( test_Hz ){
        fname = test_root+"_Hz.txt";
        outfile.open(fname.c_str());
        while( cnt < output_size ){
            outfile << z[cnt] << "\t" << cosmo->engine->get_Hz_km_s_Mpc(z[cnt]) << "\n";
            ++cnt;
        }
        outfile.close();
    }

//	testing effective EoS calculation (approximation)
    cnt = 0;
    if( test_weff ){
        double w, weff;
        fname = test_root+"_weff.txt";
        outfile.open(fname.c_str());
        while( cnt < output_size ){
            cosmo->engine->get_EoS(z[cnt],&w,&weff);
            outfile << z[cnt] << "\t" << w << "\t" << weff << "\n";
            ++cnt;
        }
        outfile.close();
    }

}
