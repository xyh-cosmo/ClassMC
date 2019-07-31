//--------------------------------------------------------------------------
//
// Description:
//     class ClassEngine :
// encapsulation of class calls
//
//
// Author List:
//    Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
//
// History (add to end):
//    creation:   ven. nov. 4 11:02:20 CET 2011
//
//  Update: I added some other useful member methods.
//-----------------------------------------------------------------------

#ifndef ClassEngine_hh
#define ClassEngine_hh

//CLASS
#include "class.h"

#include "Engine.hh"
//STD
#include <string>
#include <vector>
#include <utility>
#include <ostream>
#include <map>

#include <armadillo>

using std::string;

typedef std::map<std::string, double>  imcmc_double;    // added by Youhua Xu @ Jan-7-2017

//general utility to convert safely numerical types to string
template<typename T> std::string str(const T &x);
//specialisations
template<> std::string str (const float &x);
template<> std::string str (const double &x);
template<> std::string str (const bool &x); //"yes" or "no"
template<> std::string str (const std::string &x);

std::string str(const char* x);
//////////////////////////////////////////////////////////////////////////
//class to encapsulate CLASS parameters from any type (numerical or string)
class ClassParams {
public:

    ClassParams() {};
    ClassParams( const ClassParams& o):pars(o.pars) {};

    //use this to add a CLASS variable
    template<typename T> unsigned add(const string& key,const T& val) {
        // std::cout << "adding " << key << " with value: " << val << std::endl;
        pars.push_back(make_pair(key,str(val)));
        return pars.size();
    }

    //accesors
    inline unsigned size() const {
        return pars.size();
    }
    inline string key(const unsigned& i) const {
        return pars[i].first;
    }
    inline string value(const unsigned& i) const {
        return pars[i].second;
    }


private:
    std::vector<std::pair<string,string> > pars;
};

///////////////////////////////////////////////////////////////////////////
class ClassEngine : public Engine
{

    friend class ClassParams;

public:
    bool background_init_failed;    // default: false
    bool do_CMB;    // default true
    bool do_thermo; // default true. Needed to calculate z_rec, z_drag etc.

    //constructors
    ClassEngine(const ClassParams& pars);

    //with a class .pre file
    ClassEngine(const ClassParams& pars,const string & precision_file);

    // destructor
    ~ClassEngine();

    bool updateParValues(const std::vector<double>& par); // this version is NOT used
    bool updateParValues( std::map<std::string,std::string>& param_pair, imcmc_double& mcmc_params );

    //get value at l ( 2<l<lmax): in units = (micro-K)^2
    //don't call if FAILURE returned previously
    //throws std::execption if pb

    double getCl(Engine::cltype t,const long &l);
    void getCls(const std::vector<unsigned>& lVec, //input
                std::vector<double>& cltt,
                std::vector<double>& clte,
                std::vector<double>& clee,
                std::vector<double>& clbb);

    bool getLensing(const std::vector<unsigned>& lVec, //input
                    std::vector<double>& clphiphi,
                    std::vector<double>& cltphi,
                    std::vector<double>& clephi);

//for BAO
    inline double z_drag() const {
    //  redshift of drag epoch
        return th.z_d;
    }

    inline double rs_drag() const {
    //  comoving sound horizon at z_drag
        return th.rs_d;
    }

    inline double getTauReio() const {
        return th.tau_reio;
    }

// the following are added by Youhua Xu @ Nov-2-2015
    inline double z_rec() const { // added @Jan-19-2017
    //  redshift when recombination happens
        return th.z_rec;
    }

    inline double rs_rec() const {
    //  comoving sound horizon at z_rec
        return th.rs_rec;
    }

    inline double ra_rec() const {
    //  comoving angular diameter distance to z_dec
        return th.ra_rec;
    }

    inline double get_100theta_s() const {
        return 100.*th.rs_rec/th.ra_rec;
    }

    inline double get_a_today() const {
        return ba.a_today;
    }

    double get_Dv(double z);
    double get_Dv(double *z, double *result, int& size);  // added by YHX
    double get_Omega_m();         // added by Youhua Xu
    double get_Omega_b();         // added by Youhua Xu
    double get_Dl(double z);      // added by Youhua Xu
    double get_Da(double z);
    double get_sigma8(double z);
    double get_f(double z);
    double get_Fz(double z);
    double get_Hz(double z);
    double get_Hz_km_s_Mpc(double z);   // added by YHX
    double get_H0();                    // added by Youhua Xu
    double get_Az(double z);

    double get_rs(double z); // to get sound horizon at any redshift
    double z_drag_Hu(); // use Hu fitting formula to get z_drag and then search for the corresponding sound horizon.
    double z_rec_Hu();  // use Hu fitting formula to get CMB decoupling redshift

    void get_EoS(double z, double *w, double *weff);
    
    double get_DDE_w(double z);
	double get_DDE_weff(double z);

//  may need that
    inline int numCls() const {
        return sp.ct_size;
    };
    inline double Tcmb() const {
        return ba.T_cmb;
    }

    inline int l_max_scalars() const {
        return _lmax;
    }

//  print content of file_content
    void printFC();

//  added by Youhua Xu
    bool use_sBBN;
    double YHe_bbn;
    
//  For debug:
    void print_w0_wa_fld();

private:
    //structures class en commun
    struct file_content fc;
    struct precision pr;        /* for precision parameters */
    struct background ba;       /* for cosmological background */
    struct thermo th;           /* for thermodynamics */
    struct perturbs pt;         /* for source functions */
    struct transfers tr;        /* for transfer functions */
    struct primordial pm;       /* for primordial spectra */
    struct spectra sp;          /* for output spectra */
    struct nonlinear nl;        /* for non-linear spectra */
    struct lensing le;          /* for lensed spectra */
    struct output op;           /* for output files */

    ErrorMsg _errmsg;            /* for error messages */
    double * cl;

    //helpers
    bool dofree;
    int freeStructs();

    //call once /model
    int computeCls();

    int class_main(
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
        ErrorMsg errmsg);
    //parnames
    std::vector<std::string> parNames;

protected:


};

#endif
