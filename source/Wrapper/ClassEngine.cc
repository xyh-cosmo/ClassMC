//  NOTE: this is adapted from the original "ClassEngine.cc" in class/cpp


//--------------------------------------------------------------------------
//
// Description:
//     class ClassEngine : see header file (ClassEngine.hh) for description.
//
//------------------------------------------------------------------------
//-----------------------
// This Class's Header --
//-----------------------
#include "ClassEngine.hh"
// ------------------------
// Collaborating classes --
//-------------------------
// C++
//--------------------
#include<iostream>
#include <iomanip>
#include<string>
#include<cmath>
#include <stdexcept>
#include<sstream>
#include<numeric>
#include<cassert>

//#define DBUG

using namespace std;

template<typename T> std::string str(const T &x) {
    std::ostringstream os;
    os << x;
    return os.str();
};
//specilization
template<> std::string str (const float &x) {
    std::ostringstream os;
    os << setprecision(8) << x;
    return os.str();
}
template<> std::string str (const double &x) {
    std::ostringstream os;
    os << setprecision(16) << x;
    return os.str();
}

template<> std::string str (const bool &x) {
    {
        return x ? "yes" : "no";
    }
}

template<> std::string str (const std::string &x) {
    return x;
}

std::string str (const char* s) {
    return string(s);
}

//instanciations
template string str(const int &x);
template string str(const signed char &x);
template string str(const unsigned char &x);
template string str(const short &x);
template string str(const unsigned short &x);
template string str(const unsigned int &x);
template string str(const long &x);
template string str(const unsigned long &x);
template string str(const long long &x);
template string str(const unsigned long long &x);

//================
// Constructors //
//================
ClassEngine::ClassEngine(const ClassParams& pars): cl(0),dofree(true) {

    background_init_failed = false;
    do_CMB  = true;
    do_thermo = true;

//  prepare fp structure
    size_t n=pars.size();

    parser_init(&fc,n,"ClassMC",_errmsg);

    //config
    // =============================================================================
    // @ Dec-10-2016, note was added by Youhua Xu
    // copy ALL ClassParams from pars into file_content, this is evey important when
    // doing MCMC sampling.
    // =============================================================================
    for (size_t i=0; i<pars.size(); i++) {

        //  debug
        // cout << "p_name: " << fc.name[i] << "\t val = " << fc.value << endl;

        strcpy(fc.name[i],pars.key(i).c_str());
        strcpy(fc.value[i],pars.value(i).c_str());
        fc.read[i] = _TRUE_;    // added by XYH to avoid invalid CLASS parameter
        parNames.push_back(pars.key(i));    //store

        if (pars.key(i)=="l_max_scalars") { //identify lmax
            istringstream strstrm(pars.value(i));
            strstrm >> _lmax;
        }
    }

//     printFC();

    cout << __FILE__ << " : using lmax=" << _lmax <<endl;
    assert(_lmax>0);

//  input
    if (input_init(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_)
        throw invalid_argument(_errmsg);


//  proetction parametres mal defini
    for (size_t i=0; i<pars.size(); i++) {
        if (fc.read[i] !=_TRUE_)
            throw invalid_argument(string("(A)--invalid CLASS parameter: ")+fc.name[i]);
    }


//  calcul class
    computeCls();

    //cout <<"creating " << sp.ct_size << " arrays" <<endl;
    cl=new double[sp.ct_size];

    // printFC();
}


ClassEngine::ClassEngine(const ClassParams& pars,const string & precision_file): cl(0),dofree(true) {

    background_init_failed = false;
    do_CMB = true;
    do_thermo = true;

    struct file_content fc_precision;
    fc_precision.size = 0;
    //decode pre structure
    if (parser_read_file(const_cast<char*>(precision_file.c_str()),&fc_precision,_errmsg) == _FAILURE_) {
        throw invalid_argument(_errmsg);
    }

    //  pars
    struct file_content fc_input;
    fc_input.size = 0;
    fc_input.filename=new char[1];

    //  prepare fc par structure
    size_t n=pars.size();

    //  parser_init(&fc_input,n,"pipo",_errmsg);
    parser_init(&fc_input,n,"ClassMC",_errmsg);

    //  config
    for (size_t i=0; i<pars.size(); i++) {
        strcpy(fc_input.name[i],pars.key(i).c_str());
        strcpy(fc_input.value[i],pars.value(i).c_str());
        fc.read[i] = _TRUE_;    // added by XYH to avoid invalid CLASS parameter
        if (pars.key(i)=="l_max_scalars") {
            istringstream strstrm(pars.value(i));
            strstrm >> _lmax;
        }
    }

    cout << __FILE__ << " : using lmax=" << _lmax <<endl;
    assert(_lmax>0);

    //concatenate both
    if (parser_cat(&fc_input,&fc_precision,&fc,_errmsg) == _FAILURE_) throw invalid_argument(_errmsg);

    //parser_free(&fc_input);
    parser_free(&fc_precision);

    //input
    if (input_init(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg) == _FAILURE_)
        throw invalid_argument(_errmsg);

    //proetction parametres mal defini
    for (size_t i=0; i<pars.size(); i++) {
        if (fc.read[i] !=_TRUE_)
            throw invalid_argument(string("(B)--invalid CLASS parameter: ")+fc.name[i]);
    }

    //calcul class
    computeCls();

    //cout <<"creating " << sp.ct_size << " arrays" <<endl;
    cl=new double[sp.ct_size];

    //printFC();
}


//--------------
// Destructor --
//--------------
ClassEngine::~ClassEngine()
{
    //printFC();
    dofree && freeStructs();

    if( cl != NULL )
        delete [] cl;
}

//-----------------
// Member functions --
//-----------------
bool ClassEngine::updateParValues(const std::vector<double>& par) {

//  very time before update the cosmology, the pre-allocated memories should be freed.
    dofree && freeStructs();

    for(size_t i=0; i<par.size(); i++) {

        double val=par[i];

        //  NOTE: m_ncdm here is actually the TOTAL mass of ALL the massive neutrions
        if( strcmp(parNames[i].c_str(), "m_ncdm") == 0 ) {

            double m_ncdm_tot = par[i];
            string m_ncdm;

            if( this->ba.N_ncdm > 0 ) {

            //  TODO: improvements are needed here.
            //  we assume the massive neutrions all have the same mass
                double m_ncdm_i = m_ncdm_tot / (this->ba.N_ncdm);

                if( this->ba.N_ncdm == 1 )
                    m_ncdm = str(m_ncdm_i);
                else if( this->ba.N_ncdm == 2 )
                    m_ncdm = str(m_ncdm_i) + "," + str(m_ncdm_i);
                else if( this->ba.N_ncdm == 3 )
                    m_ncdm = str(m_ncdm_i) + "," + str(m_ncdm_i) + "," + str(m_ncdm_i);
                else
                    throw runtime_error("ClassEngine::updateParValues() ==> currently support only N_ncdm = 1,2,3");
            }

            strcpy(fc.name[i], "m_ncdm");
            strcpy(fc.value[i], m_ncdm.c_str());
        }
        else {
            strcpy(fc.name[i], parNames[i].c_str());  // this is actually not necessary ...
            strcpy(fc.value[i], str(val).c_str());
        }

        // cout << setw(15) << fc.name[i] << " = " << fc.value[i] << endl;
    }

//  ============================================================================
//  IMPORTANT NOTE:
//  all new values of the MCMC parameters are stored into the file_content
//  pointer *fc, which will be used to re-initialize all parameter settings
//  of Class, only after that true cosmological calculations will start!
//
//  All the actually parameter updates are performed inside computeCls().
//  ============================================================================

    int status = computeCls();  // real calculation is here.
    return (status==_SUCCESS_);
}


//  This new version is added by Youhua Xu @Jan-7-2017.
//  In this new version only parameters appear in mcmc_params will be updated
//  for FC
bool ClassEngine::updateParValues(  map<string,string>& param_pair,
                                    imcmc_double& mcmc_params ) {

    dofree && freeStructs();

    for( int i=0; i<fc.size; ++i ) {

        string class_par = string(fc.name[i]);

        if( param_pair.count(class_par) == 1 ) {

            string imcmc_par = param_pair[class_par];

           string old_val, new_val;
           old_val = str(fc.value[i]);

        // for our current purpose, nuetrinos' masses are fixed values, so they are not appearing here ...
           if( class_par == "A_s") {
                if( mcmc_params.count("A_s") == 1 )
                    new_val = str(mcmc_params["A_s"]);
                else if( mcmc_params.count("10^9A_s") == 1 ){
                    new_val = str(mcmc_params["10^9A_s"]*1e-9);
                }
                else if( mcmc_params.count("ln(10^10A_s)") == 1 ){
                    new_val = str(exp(mcmc_params["ln(10^10A_s)"])*1e-10);
                }
                else{ // should never happen.
                    cout << "cannot update A_s\n";
                    exit(0);
                }
            }
            else{
                new_val = str(mcmc_params[imcmc_par]);
                //  for debug
                // cout << "get new value: " << new_val << endl;
            }

            strcpy(fc.value[i], new_val.c_str());

           new_val = str(fc.value[i]);

        //  for debug
        //    if( new_val != old_val )
        //        cout << "param " << fc.name[i] << " updated from " << old_val << " to " << new_val << endl;
        }
    }


    int status = computeCls();  // real calculation is here.
    return (status==_SUCCESS_);
}

//print content of file_content
void ClassEngine::printFC() {

    printf("FILE_CONTENT SIZE=%d\n",fc.size);

    for (int i=0; i<fc.size; i++)
        printf("%d : %s = %s\n",i,fc.name[i],fc.value[i]);
}

void ClassEngine::print_w0_wa_fld(){
    //  added @ Sep-12-2017, for debugging Hyrec's (w0,wa)-settings when DDE's w(a) is used.
    cout << " >>>> DEBUG <<<<<\n";
    cout << "--> w0_fld = " << ba.w0_fld << endl;
    cout << "--> wa_fld = " << ba.wa_fld << endl;
}

int ClassEngine::class_main(
    struct file_content *pfc,
    struct precision * ppr,
    struct background * pba,
    struct thermo * pth,
    struct perturbs * ppt,
    struct transfers * ptr,
    struct primordial * ppm,
    struct spectra * psp,
    struct nonlinear * pnl,
    struct lensing * ple,
    struct output * pop,
    ErrorMsg errmsg) {

    if( use_sBBN ) {
        pth->YHe = _BBN_;
        YHe_bbn = -1;
    }

    if (input_init(pfc,ppr,pba,pth,ppt,ptr,ppm,psp,pnl,ple,pop,errmsg) == _FAILURE_) {
        printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
        dofree = false;
        return _FAILURE_;
    }


    if (background_init(ppr,pba) == _FAILURE_) {
        printf("\n\nError running background_init \n=>%s\n",pba->error_message);
        dofree = false;
        background_init_failed = true;
        return _FAILURE_;
    }

    background_init_failed = false;

    if ( do_thermo == false ){
        dofree = true;
        return _SUCCESS_;
    }

    if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
        printf("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
        background_free(&ba);
        dofree=false;
        return _FAILURE_;
    }

    if ( do_CMB == false ){ // if do_CMB is true, then do_thermo MUST also be true
        dofree = true;
        return _SUCCESS_;
    }

    YHe_bbn = pth->YHe; // make a copy of the sBBN predicted YHe.

// cout << "YES, you're calculating CMB now!!!\n";

    if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
        printf("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
        thermodynamics_free(&th);
        background_free(&ba);
        dofree=false;
        return _FAILURE_;
    }

    if (primordial_init(ppr,ppt,ppm) == _FAILURE_) {
        printf("\n\nError in primordial_init \n=>%s\n",ppm->error_message);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        dofree=false;
        return _FAILURE_;
    }

    if (nonlinear_init(ppr,pba,pth,ppt,ppm,pnl) == _FAILURE_)  {
        printf("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        dofree=false;
        return _FAILURE_;
    }

    if (transfer_init(ppr,pba,pth,ppt,pnl,ptr) == _FAILURE_) {
        printf("\n\nError in transfer_init \n=>%s\n",ptr->error_message);
        nonlinear_free(&nl);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        dofree=false;
        return _FAILURE_;
    }

    if (spectra_init(ppr,pba,ppt,ppm,pnl,ptr,psp) == _FAILURE_) {
        printf("\n\nError in spectra_init \n=>%s\n",psp->error_message);
        transfer_free(&tr);
        nonlinear_free(&nl);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        dofree=false;
        return _FAILURE_;
    }

    if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
        printf("\n\nError in lensing_init \n=>%s\n",ple->error_message);
        spectra_free(&sp);
        transfer_free(&tr);
        nonlinear_free(&nl);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        dofree=false;
        return _FAILURE_;
    }

    dofree=true;
    return _SUCCESS_;
}


int ClassEngine::computeCls() {
    return (this->class_main(&fc,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,_errmsg));
}

int ClassEngine::freeStructs() {

//  If background initialization failed, then the subsequent calculations are NOT
//  to be done, the corresponding pointers to the structures are NOT to be allocated.

    if( (do_CMB == true) && (background_init_failed == false) ) {

        if (lensing_free(&le) == _FAILURE_) {
            printf("\n\nError in spectra_free \n=>%s\n",le.error_message);
            return _FAILURE_;
        }

        if (nonlinear_free(&nl) == _FAILURE_) {
            printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
            return _FAILURE_;
        }

        if (spectra_free(&sp) == _FAILURE_) {
            printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
            return _FAILURE_;
        }

        if (primordial_free(&pm) == _FAILURE_) {
            printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
            return _FAILURE_;
        }

        if (transfer_free(&tr) == _FAILURE_) {
            printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
            return _FAILURE_;
        }

        if (perturb_free(&pt) == _FAILURE_) {
            printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
            return _FAILURE_;
        }

        if (thermodynamics_free(&th) == _FAILURE_) {
            printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
            return _FAILURE_;
        }

    }
    else if( (do_thermo == true) && (background_init_failed == false) ){
        if (thermodynamics_free(&th) == _FAILURE_) {
            printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
            return _FAILURE_;
        }
    }

    if (background_free(&ba) == _FAILURE_) {
        printf("\n\nError in background_free \n=>%s\n",ba.error_message);
        return _FAILURE_;
    }

    return _SUCCESS_;
}

double ClassEngine::getCl(Engine::cltype t,const long &l) {

    if (!dofree) throw out_of_range("no Cl available because CLASS failed");

    if (output_total_cl_at_l(&sp,&le,&op,static_cast<double>(l),cl) == _FAILURE_) {
        cerr << ">>>fail getting Cl type=" << (int)t << " @l=" << l <<endl;
        throw out_of_range(sp.error_message);
    }

    double zecl=-1;

    double tomuk=1e6*Tcmb();
    double tomuk2=tomuk*tomuk;

    switch(t)
    {
    case TT:
        (sp.has_tt==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_tt] : throw invalid_argument("no ClTT available");
        break;
    case TE:
        (sp.has_te==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_te] : throw invalid_argument("no ClTE available");
        break;
    case EE:
        (sp.has_ee==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_ee] : throw invalid_argument("no ClEE available");
        break;
    case BB:
        (sp.has_bb==_TRUE_) ? zecl=tomuk2*cl[sp.index_ct_bb] : throw invalid_argument("no ClBB available");
        break;
    case PP:
        (sp.has_pp==_TRUE_) ? zecl=cl[sp.index_ct_pp] : throw invalid_argument("no ClPhi-Phi available");
        break;
    case TP:
        (sp.has_tp==_TRUE_) ? zecl=tomuk*cl[sp.index_ct_tp] : throw invalid_argument("no ClT-Phi available");
        break;
    case EP:
        (sp.has_ep==_TRUE_) ? zecl=tomuk*cl[sp.index_ct_ep] : throw invalid_argument("no ClE-Phi available");
        break;
    }

    return zecl;

}

void ClassEngine::getCls(const std::vector<unsigned>& lvec, //input
                         std::vector<double>& cltt,
                         std::vector<double>& clte,
                         std::vector<double>& clee,
                         std::vector<double>& clbb)
{
    cltt.resize(lvec.size());
    clte.resize(lvec.size());
    clee.resize(lvec.size());
    clbb.resize(lvec.size());

    for (size_t i=0; i<lvec.size(); i++) {
        try {
            cltt[i]=getCl(ClassEngine::TT,lvec[i]);
            clte[i]=getCl(ClassEngine::TE,lvec[i]);
            clee[i]=getCl(ClassEngine::EE,lvec[i]);
            clbb[i]=getCl(ClassEngine::BB,lvec[i]);
        }
        catch(exception &e) {
            throw e;
        }
    }

}

bool ClassEngine::getLensing(const std::vector<unsigned>& lvec, //input
                             std::vector<double>& clpp,
                             std::vector<double>& cltp,
                             std::vector<double>& clep  ) {

    clpp.resize(lvec.size());
    cltp.resize(lvec.size());
    clep.resize(lvec.size());

    for (size_t i=0; i<lvec.size(); i++) {
        try {
            clpp[i]=getCl(ClassEngine::PP,lvec[i]);
            cltp[i]=getCl(ClassEngine::TP,lvec[i]);
            clep[i]=getCl(ClassEngine::EP,lvec[i]);
        }
        catch(exception &e) {
            cout << "plantage!" << endl;
            cout << __FILE__ << e.what() << endl;
            return false;
        }
    }

    return true;
}


double ClassEngine::get_f(double z) {
    double tau;
    int index;
    double *pvecback;   //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);    //call to fill pvecback
    double f = pvecback[ba.index_bg_f];
    free(pvecback);
    return f;
}


double ClassEngine::get_sigma8(double z) {

    int index;
    double tau;
    double *pvecback;
    double sigma8 = 0.;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback

    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
    //background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback);
    spectra_sigma(&ba,&pm,&sp,8./ba.h,z,&sigma8);

//    cout << "@ tau = " << setprecision(12) << setw(12) << tau << endl;

    // debug
/*
    z = 0.;
    while( z < 1.0 ){
        spectra_sigma(&ba,&pm,&sp,8./ba.h,z,&sigma8);
        cout << "z = " << z << ", sigma8 = " << sigma8 << endl;
        z += 0.1;
    }
    exit(0);
*/

    free(pvecback);
    return sigma8;
}

// ATTENTION FONCTION BIDON - GET omegam ! ------------------
double ClassEngine::get_Az(double z) {
    double Dv = get_Dv(z);
    // A(z)=100DV(z)sqrt(~mh2)/cz
    double omega_bidon = 0.12 ;
    // double Az = 100.*Dv*sqrt(omega_bidon)/(3.e8*z); // is there speed of light somewhere ?
    double Az = 100.*Dv*sqrt(omega_bidon)/(_c_*z); // _c_ is defined in background.h
    return Az;
}

void ClassEngine::get_EoS(double z, double *w, double *weff){
    background_DDE_get_EoS(&ba,z,w,weff);
}

double ClassEngine::get_Dv(double z) {

    int index;
    double tau;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

    double H_z=pvecback[ba.index_bg_H];
    double D_ang=pvecback[ba.index_bg_ang_distance];
    double D_v;
    D_v=pow(D_ang*(1+z),2)*z/H_z;
    D_v=pow(D_v,1./3.);

    free(pvecback);

    return D_v;
}

double ClassEngine::get_Fz(double z) {

    int index;
    double tau;
    double *pvecback;

    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);

    double H_z=pvecback[ba.index_bg_H];
    double D_ang=pvecback[ba.index_bg_ang_distance];
    // double F_z = (1.+z) * D_ang * H_z /(3.e8) ; // is there speed of light somewhere ?
    double F_z = (1.+z) * D_ang * H_z /_c_ ;
    free(pvecback);

    return F_z;
}

double ClassEngine::get_Hz(double z) {  // this default version use unit 1/Mpc

    int index;
    double tau;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
    double Hz = pvecback[ba.index_bg_H];
    free(pvecback);
    return Hz;
}

double ClassEngine::get_Hz_km_s_Mpc(double z) {

    int index;
    double tau;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
    double Hz = pvecback[ba.index_bg_H]*_c_*1e-3;
    free(pvecback);
    return Hz;
}

double ClassEngine::get_H0() {
    return ba.H0*_c_*1e-3;
}

//  angular diameter distance
double ClassEngine::get_Da(double z) {
    double tau;
    int index;
    double *pvecback;
    //transform redshift in conformal time
    background_tau_of_z(&ba,z,&tau);

    //pvecback must be allocated
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));

    //call to fill pvecback
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
    double D_ang=pvecback[ba.index_bg_ang_distance];
    free(pvecback);

    return D_ang;
}

//  luminosity distance
double ClassEngine::get_Dl(double z) {
    int index;
    double tau, *pvecback;
    background_tau_of_z(&ba,z,&tau);
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
    double Dl = pvecback[ba.index_bg_lum_distance];
    free(pvecback);
    return Dl;
}

double ClassEngine::get_Omega_m() {
// NOTE: I define Omega0_m = Omega0_b + Omega0_c
    return ba.Omega0_b + ba.Omega0_cdm;
}

double ClassEngine::get_Omega_b() {
    // cout << "Omega0_b = " << ba.Omega0_b << endl;
    // exit(0);
    return ba.Omega0_b;
}


double ClassEngine::z_drag_Hu(){
    double h = 0.01 * get_H0();
    double omegab = ba.Omega0_b*h*h;
    double omegam = ba.Omega0_m*h*h;
    // double omegam = get_Omega_m()*h*h;
    double b1 = 0.313*pow(omegam, -0.419)*(1+0.607*pow(omegam, 0.674));
    double b2 = 0.238*pow(omegam, 0.223);
    double zdrag = 1291*pow(omegam,0.251)*(1+b1*pow(omegab,b2))/(1+0.659*pow(omegam,0.828));
    return zdrag;
}


double ClassEngine::z_rec_Hu(){
    double h = 0.01 * get_H0();
    double omegab = ba.Omega0_b*h*h;
    double omegam = ba.Omega0_m*h*h;
    double g1 = 0.0783*pow(omegab, -0.238)/(1+39.5*pow(omegab, 0.763));
    double g2 = 0.560/(1+21.1*pow(omegab, 1.81));
    double zrec = 1048.*(1+0.00124*pow(omegab,-0.738))*(1+g1*pow(omegam,g2));
	return zrec;
}


double ClassEngine::get_rs(double z){
    // printf("class_z_drag = %g\n",th.z_d);
    // printf("Hu_z_drag    = %g\n",z);
    int index;
    double tau, *pvecback;
    background_tau_of_z(&ba,z,&tau);
    pvecback=(double *)malloc(ba.bg_size*sizeof(double));
    background_at_tau(&ba,tau,ba.long_info,ba.inter_normal, &index, pvecback);
    double rs = pvecback[ba.index_bg_rs];
    free(pvecback);
    return rs;
}

double ClassEngine::get_DDE_w(double z){
	double w,weff;
	background_DDE_get_EoS( &ba,z,&w,&weff);
	return w;
}

double ClassEngine::get_DDE_weff(double z){
	double w,weff;
	background_DDE_get_EoS( &ba,z,&w,&weff);
	return weff;
}
