#include "ClassMC.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

// void get_plank_extra_params( CosmoTheory*           model,
//                              Data_Planck2015*       data,
//                              imcmc_vector_string&   params ){
//     Data_Planck2015*    data_plk = static_cast<Data_Planck2015*>(data);

//     for( size_t i=0; i<data_plk->extra_names.size(); ++i ){
//         params.push_back(data_plk->extra_names[i]);
//     }
// }

Data_Planck2015::Data_Planck2015(){
    lmax_lens[0] = -1;

    for( int i=0; i<6; ++i ){
        lmax_highl[i]   = -1;
        has_cl_highl[i] = 0;

        lmax_lowP[i]    = -1;
        has_cl_lowP[i]  = 0;

        lmax_lens[i+1] = -1;
    }
}

bool Data_Planck2015::Init( string& dataset ){

    data_info.GetInfo(dataset);

    bool init_status = false;
    _err    = initError();
    err     = &_err;

    sprintf(clnames[0],"TT");
    sprintf(clnames[1],"EE");
    sprintf(clnames[2],"BB");
    sprintf(clnames[3],"TE");
    sprintf(clnames[4],"TB");    //    assumed to be zero
    sprintf(clnames[5],"EB");    //    assumed to be zero

    string clikfile_highl, clikfile_lowP, clikfile_lens;

    if( Read::Has_Key_in_File(dataset, "clikfile_highl") )
        clikfile_highl = Read::Read_String_from_File(dataset, "clikfile_highl");
    else
        clikfile_highl = "";

    if( Read::Has_Key_in_File(dataset, "clikfile_lowP") )
        clikfile_lowP = Read::Read_String_from_File(dataset, "clikfile_lowP");
    else
        clikfile_lowP = "";

    if( Read::Has_Key_in_File(dataset, "clikfile_lens") )
        clikfile_lens = Read::Read_String_from_File(dataset, "clikfile_lens");
    else
        clikfile_lens = "";

    has_highl   = (clikfile_highl != "");
    has_lowP    = (clikfile_lowP != "");
    has_lens    = (clikfile_lens != "");

    char clikfile_highl_[clikfile_highl.size()+1];
    char clikfile_lowP_[clikfile_lowP.size()+1];
    char clikfile_lens_[clikfile_lens.size()+1];

    strcpy(clikfile_highl_, clikfile_highl.c_str());
    strcpy(clikfile_lowP_, clikfile_lowP.c_str());
    strcpy(clikfile_lens_, clikfile_lens.c_str());

    extra_names.clear();
    bool same_name = false;

    if( has_highl ){//    init highl

        highl    = clik_init( clikfile_highl_, err);
        quitOnError(*err,__LINE__,stderr);

        clik_get_has_cl( highl, has_cl_highl, err );
        quitOnError(*err,__LINE__,stderr);

        clik_get_lmax( highl, lmax_highl, err );
        quitOnError(*err,__LINE__,stderr);

        // fprintf(stdout,"Likelihood use Cl\n");
        // for(cli=0; cli<6; cli++) {
        //     if( has_cl_highl[cli]==1 ){
        //         fprintf(stdout,"  %s from l=0 to l=%d (incl.)\n", clnames[cli], lmax_highl[cli]);
        //     }
        // }

        nextra_highl = clik_get_extra_parameter_names( highl, &names_highl, err);
        quitOnError(*err,__LINE__,stderr);
        if ( nextra_highl > 0 ) {
            // fprintf(stdout,"With %d extra parameters\n", nextra_highl);
            // for( i=0; i<nextra_highl; i++)
            //       fprintf(stdout,"  %s\n", names_highl[i]);
            //
            ndim_highl = nextra_highl;
        }
        else
            ndim_highl = 0;

        for( cli=0; cli<6; cli++ ) {
            ndim_highl += (lmax_highl[cli] + 1);
        }

        if( ndim_highl > 0 ){
            cl_and_pars_highl = new double[ndim_highl];

            for( int l=0; l<ndim_highl; ++l )
                cl_and_pars_highl[l] = 0.0;
        }
        else
            throw runtime_error(" cl_and_pars_highl has zero size ");

    //  add names_highl into extra_names
        if( nextra_highl > 0 ){

            for( int n=0; n<nextra_highl; ++n ){
                string temp_name(names_highl[n]);

                same_name = false;

                if( extra_names.size() == 0 ){
                    extra_names.push_back(temp_name);
                }
                else{
                    for(size_t k=0; k<extra_names.size(); ++k){
                        if( Read::SameStrings(temp_name, extra_names[k]) ){
                            same_name = true;
                            break;
                        }
                    }
                    if( same_name == false )
                        extra_names.push_back(temp_name);
                    else
                        sprintf(names_highl[n], "%s", temp_name.c_str());   // make sure all likelihoods use exactly the same parameter names
                }

                // check if there is A_Planck
                if( Read::SameStrings(temp_name, "A_Planck") ){
					has_A_Planck = true;
				}
            }
        }

    }

    // init lowP
    if( has_lowP ){

        lowP    = clik_init( clikfile_lowP_, err);
        quitOnError(*err,__LINE__,stderr);

        clik_get_has_cl( lowP, has_cl_lowP, err );
        quitOnError(*err,__LINE__,stderr);

        clik_get_lmax( lowP, lmax_lowP, err );
        quitOnError(*err,__LINE__,stderr);

        // fprintf(stdout,"Likelihood use Cl\n");
        // for(cli=0; cli<6; cli++) {
        //     if( has_cl_lowP[cli]==1 ){
        //         fprintf(stdout,"  %s from l=0 to l=%d (incl.)\n", clnames[cli], lmax_lowP[cli]);
        //     }
        // }

        nextra_lowP = clik_get_extra_parameter_names( lowP, &names_lowP, err);
        quitOnError(*err,__LINE__,stderr);

        if ( nextra_lowP > 0 ) {
            // fprintf(stdout,"With %d extra parameters\n", nextra_lowP);
            // for( i=0; i<nextra_lowP; i++)
            //       fprintf(stdout,"  %s\n", names_lowP[i]);

            ndim_lowP = nextra_lowP;
        }
        else
            ndim_lowP = 0;

        for( cli=0; cli<6; cli++ ) {
            ndim_lowP += (lmax_lowP[cli] + 1);
        }

        if( ndim_lowP > 0 ){
            cl_and_pars_lowP = new double[ndim_lowP];

            for( int l=0; l<ndim_lowP; ++l )
                cl_and_pars_lowP[l] = 0.0;
        }
        else
            throw runtime_error(" cl_and_pars_lowP has zero size ");

        if( nextra_lowP > 0 ){
            for( int n=0; n<nextra_lowP; ++n ){

                string temp_name(names_lowP[n]);

                same_name = false;
                if( extra_names.size() == 0 ){
                    extra_names.push_back(temp_name);
                }
                else{
                    for(size_t k=0; k<extra_names.size(); ++k){
                        if( Read::SameStrings(temp_name, extra_names[k]) ){
                            same_name = true;
                            break;
                        }
                    }

                    if( same_name == false )
                        extra_names.push_back(temp_name);
                    else
                        sprintf(names_lowP[n], "%s", temp_name.c_str());   // make sure all likelihoods use exactly the same parameter names
                }

                // check if there is A_Planck
                if( Read::SameStrings(temp_name, "A_Planck") ){
					has_A_Planck = true;
					if( temp_name == "A_planck" ){
                    //  planck commander likelihood use 'A_planck' instead of 'A_Planck' ... well, who wrote commander???
						sprintf(names_lowP[n], "%s", "A_Planck");
					}
				}
            }
        }
    }

    if( has_lens ){//    init lens

        cout << "#---> initializing lensing likelihood <---\n";

        lens     = clik_lensing_init( clikfile_lens_, err);
        quitOnError(*err,__LINE__,stderr);

        cout << "#---> geting lmaxs[] of lening likelihood <---\n";
        clik_lensing_get_lmaxs( lens, lmax_lens, err );
        quitOnError(*err,__LINE__,stderr);

        fprintf(stdout,"Likelihood use Cl\n");
        fprintf(stdout,"  PP from l=0 to l=%d (incl.)\n", lmax_lens[0]);

        fprintf(stdout,"Likelihood use Cl\n");
        for(cli=1; cli<7; cli++) {    //    be careful here!!!
            if( lmax_lens[cli] > 0 ){
                fprintf(stdout,"  %s from l=0 to l=%d (incl.)\n", clnames[cli-1], lmax_lens[cli]);
            }
        }

        nextra_lens = clik_lensing_get_extra_parameter_names( lens, &names_lens, err);
        quitOnError(*err,__LINE__,stderr);

        if ( nextra_lens > 0 ) {
            fprintf(stdout,"#  Planck2015 lensing likelihood with %d extra parameters\n", nextra_lens);
            for( i=0; i<nextra_lens; i++)
                  fprintf(stdout,"  %s\n", names_lens[i]);

            ndim_lens = nextra_lens;    //    however, nextra_lens = -1
        }
        else
            ndim_lens = 0;

        for( cli=0; cli<7; cli++ ) {    //    BE CAREFULL !!!
            ndim_lens += (lmax_lens[cli] + 1);
        }

        if( ndim_lens > 0 ){
            cl_and_pars_lens = new double[ndim_lens];

            for( int l=0; l<ndim_lens; ++l )
                cl_and_pars_lens[l] = 0.0;
        }
        else
            throw runtime_error(" cl_and_pars has zero size ");

        if( nextra_lens > 0 ){
            for( int n=0; n<nextra_lens; ++n ){
                string temp_name(names_lens[n]);

                same_name = false;
                if( extra_names.size() == 0 )
                    extra_names.push_back(temp_name);
                else{
                    for(size_t k=0; k<extra_names.size(); ++k){
                        if( Read::SameStrings(temp_name, extra_names[k]) ){
                            same_name = true;
                            break;
                        }
                    }

                    if( same_name == false )
                        extra_names.push_back(temp_name);
                    else
                        sprintf(names_lens[n], "%s", temp_name.c_str());
                }

                // check if there is A_Planck
                if( Read::SameStrings(temp_name, "A_Planck") ){
					has_A_Planck = true;
					sprintf(names_lens[n], "%s", "A_Planck");   // make sure all likelihoods use exactly the same parameter names
				}
            }
        }
    }

    //  ============
    //  A_Planck
    //  ============

    if( Read::Has_Key_in_File(dataset, "A_Planck_mean") )
        A_Planck_mean   = Read::Read_Double_from_File(dataset, "A_Planck_mean");
    else
        A_Planck_mean   = 1.0;      //  default value

    if( Read::Has_Key_in_File(dataset, "A_Planck_std") )
        A_Planck_std    = Read::Read_Double_from_File(dataset, "A_Planck_std");
    else
        A_Planck_std    = 0.0025;    //  default value

//  ================================================================================
//  add extra parameters needed by clik likelihoods. These parameters will be laster
//  added to the MCMC_Params inside datalist.
    for( size_t i=0; i<extra_names.size(); ++i ){
        data_info.Add_MCMC_Param(extra_names[i]);
        cout << "adding clik extra parameter: " << extra_names[i] << endl;
    }

    init_status = has_highl || has_lowP || has_lens;
    return init_status;
}

Data_Planck2015::~Data_Planck2015(){

    if( has_highl ){
        delete[] cl_and_pars_highl;
        clik_cleanup(&highl);

        if( nextra_highl > 0 )
            free(names_highl);
    }

    if( has_lowP ){
        delete[] cl_and_pars_lowP;
        clik_cleanup(&lowP);

        if( nextra_lowP > 0 )
            free(names_lowP);
    }

    if( has_lens ){
        delete[] cl_and_pars_lens;
        clik_lensing_cleanup(&lens);

        if( nextra_lens > 0 )
            free(names_lens);
    }
}

//========================================================
//enum cltype {TT=0,EE,TE,BB,PP,TP,EP}, defined in class
//    TT  = 0
//    EE  = 1
//    BB  = 3
//    TE  = 2
//========================================================

void Data_Planck2015::get_Cl_tasks( bool& has_tCl, bool& has_pCl, bool& has_lCl ){

    has_tCl = has_pCl = has_lCl = false;

    //  tCl
    {
        if( has_cl_highl[0] == 1 )
            has_tCl = true;

        if( has_cl_highl[1] || has_cl_highl[3] || has_cl_highl[3] || has_cl_highl[4] || has_cl_highl[5] )
            has_pCl = true;
    }

    //  pCl
    {
        if( has_cl_lowP[0] == 1 )
            has_tCl = true;

        if( has_cl_lowP[1] || has_cl_lowP[3] || has_cl_lowP[3] || has_cl_lowP[4] || has_cl_lowP[5] )
            has_pCl = true;
    }

    //  lCl
    {
        if( lmax_lens[0] > 0 )
            has_lCl = true;

        if( lmax_lens[1] > 0 )
            has_tCl = true;

        if( lmax_lens[2]>0 || lmax_lens[3]>0 || lmax_lens[4]>0 || lmax_lens[5]>0 || lmax_lens[6]>0 )
            has_pCl = true;
    }
}

void Data_Planck2015::get_cls( CosmoTheory* theory ){
    if( has_highl ){
        // cout << "getting highl Cls ....\n";
        get_cls_highl( theory );
    }

    if( has_lowP ) {
        // cout << "getting lowP Cls ...\n";
        get_cls_lowP( theory );
    }

    if( has_lens ) {
        // cout << "getting lens Cls ...\n";
        get_cls_lens( theory );
    }
}

void Data_Planck2015::get_cls_highl( CosmoTheory* theory ){

    nall_highl    = 0;

    for( int l=0; l<=lmax_highl[0]; ++l ){//    get TT from class
        if( l >= 2 )
            cl_and_pars_highl[nall_highl] = theory->engine->getCl(theory->engine->TT, l);
        else
            cl_and_pars_highl[nall_highl] = 0;

        ++nall_highl;
    }

    for( int l=0; l<=lmax_highl[1]; ++l ){//    get EE from class
        if( l >= 2 )
            cl_and_pars_highl[nall_highl] = theory->engine->getCl(theory->engine->EE, l);
        else
            cl_and_pars_highl[nall_highl] = 0;

        ++nall_highl;
    }

    for( int l=0; l<=lmax_highl[2]; ++l ){//    get BB from class
        if( l >= 2 )
            cl_and_pars_highl[nall_highl] = theory->engine->getCl(theory->engine->BB, l);
        else
            cl_and_pars_highl[nall_highl] = 0;

        ++nall_highl;
    }

    for( int l=0; l<=lmax_highl[3]; ++l ){//    get TE from class
        if( l >= 2 )
            cl_and_pars_highl[nall_highl] = theory->engine->getCl(theory->engine->TE, l);
        else
            cl_and_pars_highl[nall_highl] = 0;

        ++nall_highl;
    }

#ifdef __CLASSMC_DEBUG__
    if( (nall_highl + nextra_highl) != ndim_highl ){
        cout << "**** nall_highl   = " << nall_highl << endl;
        cout << "**** nextra_highl = " << nextra_highl << endl;
        cout << "**** ndim_highl   = " << ndim_highl << endl;
        throw runtime_error("**** error when filling cl_and_pars_highl ****");
    }
#endif
}

void Data_Planck2015::get_cls_lowP( CosmoTheory* theory ){

    nall_lowP    = 0;

    for( int l=0; l<=lmax_lowP[0]; ++l ){//    get TT from class
        if( l >= 2 )
            cl_and_pars_lowP[nall_lowP] = theory->engine->getCl(theory->engine->TT, l);
        else
            cl_and_pars_lowP[nall_lowP] = 0;

        ++nall_lowP;
    }

    for( int l=0; l<=lmax_lowP[1]; ++l ){//    get EE from class
        if( l >= 2 )
            cl_and_pars_lowP[nall_lowP]    = theory->engine->getCl(theory->engine->EE, l);
        else
            cl_and_pars_lowP[nall_lowP] = 0;

        ++nall_lowP;
    }

    for( int l=0; l<=lmax_lowP[2]; ++l ){//    get BB from class
        if( l >= 2 )
            cl_and_pars_lowP[nall_lowP] = theory->engine->getCl(theory->engine->BB, l);
        else
            cl_and_pars_lowP[nall_lowP] = 0;

        ++nall_lowP;
    }

    for( int l=0; l<=lmax_lowP[3]; ++l ){//    get TE from class
        if( l >= 2 )
            cl_and_pars_lowP[nall_lowP] = theory->engine->getCl(theory->engine->TE, l);
        else
            cl_and_pars_lowP[nall_lowP] = 0;

        ++nall_lowP;
    }

#ifdef __CLASSMC_DEBUG__
    if( (nall_lowP + nextra_lowP) != ndim_lowP ){
        cout << "**** nall_lowP   = " << nall_lowP << endl;
        cout << "**** nextra_lowP = " << nextra_lowP << endl;
        cout << "**** ndim_lowP   = " << ndim_lowP << endl;
        throw runtime_error("**** error when filling cl_and_pars_lowP ****");
    }
#endif
}

void Data_Planck2015::get_cls_lens( CosmoTheory* theory ){

    //=================================
    //    currently only PP is included
    //=================================

    nall_lens    = 0;

    for( int l=0; l<=lmax_lens[0]; ++l ){    //    get PhiPhi
        if( l >= 2 )
            cl_and_pars_lens[nall_lens] = theory->engine->getCl(theory->engine->PP, l);
        else
            cl_and_pars_lens[nall_lens] = 0;

        ++nall_lens;
    }

    for( int l=0; l<=lmax_lens[1]; ++l ){    //    get TT from class
        if( l >= 2 )
            cl_and_pars_lens[nall_lens] = theory->engine->getCl(theory->engine->TT, l);
        else
            cl_and_pars_lens[nall_lens] = 0;

        ++nall_lens;
    }

    for( int l=0; l<=lmax_lens[2]; ++l ){    //    get EE from class
        if( l >= 2 )
            cl_and_pars_lens[nall_lens]    = theory->engine->getCl(theory->engine->EE, l);
        else
            cl_and_pars_lens[nall_lens] = 0;

        ++nall_lens;
    }

    for( int l=0; l<=lmax_lens[3]; ++l ){    //    get BB from class
        if( l >= 2 )
            cl_and_pars_lens[nall_lens] = theory->engine->getCl(theory->engine->BB, l);
        else
            cl_and_pars_lens[nall_lens] = 0;

        ++nall_lens;
    }

    for( int l=0; l<=lmax_lens[4]; ++l ){    //    get TE from class
        if( l >= 2 )
            cl_and_pars_lens[nall_lens] = theory->engine->getCl(theory->engine->TE, l);
        else
            cl_and_pars_lens[nall_lens] = 0;
        
        ++nall_lens;
    }


#ifdef __CLASSMC_DEBUG__
    if( (nall_lens + nextra_lens) != ndim_lens ){
//    currently lensing likelihood has no extra parameter, don't believe clik_print
        cout << "**** nall_lens   = " << nall_lens << endl;
        cout << "**** nextra_lens = " << nextra_lens << endl;
        cout << "**** ndim_lens   = " << ndim_lens << endl;
        throw runtime_error("**** error when filling cl_and_pars_lens ****");
    }
#endif

}

//  ================================
//  MUST be called after get_cls()
//  ================================
void Data_Planck2015::update_extra_params( imcmc_double& full_param ){

    if( has_highl && (nextra_highl>0) ){
        for( int i=0; i<nextra_highl; ++i ){
            cl_and_pars_highl[nall_highl] = full_param[names_highl[i]];
            ++nall_highl;
        }
    }

    if( has_lowP && (nextra_lowP>0) ){
        for( int i=0; i<nextra_lowP; ++i ){
            cl_and_pars_lowP[nall_lowP] = full_param[names_lowP[i]];
            ++nall_lowP;
        }
    }

    if( has_lens && (nextra_lens>0) ){
        for( int i=0; i<nextra_lens; ++i ){
            cl_and_pars_lens[nall_lens] = full_param[names_lens[i]];
            ++nall_lens;
        }
    }
    
//  compute the Planck recommended prior on A_Planck.
    double A_Planck = -1.0;
    loglike_A_Planck = 0.0;

    if( has_A_Planck ){

        if( full_param.count("A_planck") == 1){
            A_Planck = full_param["A_planck"];
        }
        else if( full_param.count("A_Planck") == 1 ){
            A_Planck = full_param["A_Planck"];
        }

        if( A_Planck > 0 ){ // make sure got a sensitive value for A_Planck
            loglike_A_Planck = -0.5 * pow( (A_Planck - A_Planck_mean)/A_Planck_std, 2.);
           // cout << "debug loglike_A_Planck = " << loglike_A_Planck << endl;
        }
    }
}

void Data_Planck2015::compute_chisq(){

    // lens likelihood is problematic ..

    loglike_highl = loglike_lowP = loglike_lens = 0.0;

    if( has_highl )
        loglike_highl = clik_compute( highl, cl_and_pars_highl, err );

    if( has_lowP )
        loglike_lowP  = clik_compute( lowP, cl_and_pars_lowP, err );

    if( has_lens )
        loglike_lens  = clik_lensing_compute( lens, cl_and_pars_lens, err );

    cout << "loglike_highl = " << loglike_highl << endl;
    cout << "loglike_lowP  = " << loglike_lowP << endl;
    cout << "loglike_lens  = " << loglike_lens << endl;

    loglike_total = (loglike_highl + loglike_lowP + loglike_lens + loglike_A_Planck);
//    exit(0);
}
