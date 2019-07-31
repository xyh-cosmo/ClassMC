#include "ClassMC.hpp"
#include "WMAP7.hpp"
#include <sstream>

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;


double Likelihood_WMAP7(imcmc_double	&full_params,
                        double			&lndet,
                        double			&chisq,
                        void			*model,
                        void			*data,
                        istate			&state ) {

    lndet = chisq = 0.0;

    CosmoTheory *cosmo = static_cast<CosmoTheory*>(model);
    Data_WMAP7	*WMAP7 = static_cast<Data_WMAP7*>(data);

//	get Cls
    double *cltt, *clte, *clee, *clbb;
#if defined(USE_LOWELL_TBEB) || defined(USE_HIGHELL_TB)
    double *cltb, *cleb;
#endif

    int size=WMAP7->lmax-2+1;

    cltt = new double[size];
    clte = new double[size];
    clee = new double[size];
    clbb = new double[size];

#if defined(USE_LOWELL_TBEB) || defined(USE_HIGHELL_TB)
    cltb = new double[size];
    cleb = new double[size];
#endif

    for( int l=2; l<=WMAP7->lmax; ++l ) {

        double factor = 0.5*l*(l+1.0)/_PI_;

        cltt[l-2] = cosmo->engine->getCl(cosmo->engine->TT, l)*factor;
        clte[l-2] = cosmo->engine->getCl(cosmo->engine->TE, l)*factor;
        clee[l-2] = cosmo->engine->getCl(cosmo->engine->EE, l)*factor;
        clbb[l-2] = cosmo->engine->getCl(cosmo->engine->BB, l)*factor;
#if defined(USE_LOWELL_TBEB)// || defined(USE_HIGHELL_TB)
        cltb[l-2] = cosmo->engine->getCl(cosmo->engine->TB, l)*factor;
        cleb[l-2] = cosmo->engine->getCl(cosmo->engine->EB, l)*factor;
#endif
    }

//	compute WMAP7 likelihood
#ifdef USE_HIGHELL_TB
    int num_WMAP = 10; // number of individual chi2 terms in likelihood
#else
    int num_WMAP = 8; // number of individual chi2 terms in likelihood
#endif

    double *like = new double[num_WMAP];


// NOTE: intel fortran compiler adds extra under line '_' to the fun-name-end, but
//		gfortran adds '__' to the fun-name-head
#ifdef USE_INTEL
    #if defined(USE_LOWELL_TBEB) || defined(USE_HIGHELL_TB)
        wmap_likelihood_7yr_mp_wmap_likelihood_compute_(cltt,clte,cltb,clee,cleb,clbb,like);
    #else
        wmap_likelihood_7yr_mp_wmap_likelihood_compute_(cltt,clte,clee,clbb,like);
    #endif
#elif USE_GFORTRAN
    #if defined(USE_LOWELL_TBEB) || defined(USE_HIGHELL_TB)
        __wmap_likelihood_7yr_MOD_wmap_likelihood_compute(cltt,clte,cltb,clee,cleb,clbb,like);
    #else
        __wmap_likelihood_7yr_MOD_wmap_likelihood_compute(cltt,clte,clee,clbb,like);
    #endif
#endif

    double like_tot = 0.0;
    int like_idx=0;
    while ( like_idx < num_WMAP ) {
        like_tot += like[like_idx];
        ++like_idx;
    }

//  ======================
//  free allocated memory
    delete[] like;
    delete[] cltt;
    delete[] clte;
    delete[] clee;
    delete[] clbb;

#if defined(USE_LOWELL_TBEB)// || defined(USE_HIGHELL_TB)
    delete[] cltb;
    delete[] cleb;
#endif

    chisq += 2.0*like_tot;
    return lndet - 0.5*chisq;
}
