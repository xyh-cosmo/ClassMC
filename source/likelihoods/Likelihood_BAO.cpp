
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "CosmoTheory.hpp"
#include "BAO.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;

void Compute_BAO_Vals(	CosmoTheory* cosmo,
                        double *bao_z,
                        double *bao_prediction,
                        int& bao_size,
                        int& bao_type,
                        bool& use_Hu_fitting  ) {
    // time_t _time_;
    // DEBUG_TIME_START(_time_);

    double rs;
    double Da,Dv;
    double H0,Om;

    switch(bao_type) {
        case _bao_rs_over_DV_:
            for( int i=0; i<bao_size; ++i ) {
                Dv = cosmo->engine->get_Dv(bao_z[i]);
                use_Hu_fitting ? \
                    rs = cosmo->engine->get_rs(cosmo->engine->z_drag_Hu()) : \
                    rs = cosmo->engine->rs_drag();
                bao_prediction[i] = rs/Dv;
            }
            break;
        case _bao_A_:
            for( int i=0; i<bao_size; ++i ) {
                Dv = cosmo->engine->get_Dv(bao_z[i]);		// unit = Mpc
                Om = cosmo->engine->get_Omega_m();
                H0 = cosmo->engine->get_H0()/(1e-3*_c_); // convert to 1/Mpc
                bao_prediction[i] = H0*Dv*sqrt(Om)/bao_z[i];

                // std::cout << "## z = " << bao_z[i] << ", Dv = " << Dv << ", Om = " << Om << ", H0 = " << H0 << std::endl;
            }
            break;
        case _bao_DV_over_rs_:
            for( int i=0; i<bao_size; ++i ) {
                Dv = cosmo->engine->get_Dv(bao_z[i]);
                use_Hu_fitting ? \
                    rs = cosmo->engine->get_rs(cosmo->engine->z_drag_Hu()) : \
                    rs = cosmo->engine->rs_drag();
                bao_prediction[i] = Dv/rs;

                // std::cout << "## z = " << bao_z[i] << ", Dv = " << Dv << ", rs = " << rs << std::endl;
            }
            break;
        case _bao_DV_:
            for( int i=0; i<bao_size; ++i ) {
                bao_prediction[i] = cosmo->engine->get_Dv(bao_z[i]);
                // std::cout << "## z = " << bao_z[i] << ", Dv = " << bao_prediction[i] << std::endl;
            }
            break;
        case _bao_DA_over_rs_:
            for( int i=0; i<bao_size; ++i ) {
                Da = cosmo->engine->get_Da(bao_z[i]);
                use_Hu_fitting ? \
                    rs = cosmo->engine->get_rs(cosmo->engine->z_drag_Hu()) : \
                    rs = cosmo->engine->rs_drag();
                bao_prediction[i] = Da/rs;

                // std::cout << "## z = " << bao_z[i] << ", Da = " << Da << ", rs = " << rs << std::endl;
            }
            break;
        case _bao_DA_:
            for( int i=0; i<bao_size; ++i ) {
                Da = cosmo->engine->get_Da(bao_z[i]);
                use_Hu_fitting ? \
                    rs = cosmo->engine->get_rs(cosmo->engine->z_drag_Hu()) : \
                    rs = cosmo->engine->rs_drag();
                bao_prediction[i] = Da;

                // std::cout << "## z = " << bao_z[i] << ", Da = " << Da << ", rs = " << rs << std::endl;
            }
            break;
        case _bao_DR12_zhao_: // the format is adopted to use DR12 Zhao 2016 BAO results
            int nz = bao_size/2;
            for( int i=0; i<nz; ++i ) {
                double Hz = cosmo->engine->get_Hz_km_s_Mpc(bao_z[i]);
                Da = cosmo->engine->get_Da(bao_z[i]);
                use_Hu_fitting ? \
                    rs = cosmo->engine->get_rs(cosmo->engine->z_drag_Hu()) : \
                    rs = cosmo->engine->rs_drag();
                bao_prediction[i] = Da/rs;
				bao_prediction[nz+i] = Hz*rs;
//                 std::cout << "## z = " << bao_z[i] << ", Hz = " << Hz << ", rs = " << rs << std::endl;
            }
            break;
    }

    // DEBUG_TIME_END(_time_);
}

double Compute_BAO_Chi2(double *bao_data,
                        double *bao_prediction,
                        double *bao_err,
                        double *bao_icov,
                        int& n_bao) {

    // time_t _time_;
    // DEBUG_TIME_START(_time_);

    double chi2 = 0.0;

    if( n_bao == 1 ) {
        chi2 = pow( (bao_data[0]-bao_prediction[0])/bao_err[0] ,2);
    }
    else if( n_bao > 1 ) {
        for( int i=0; i<n_bao; ++i ) {
            for( int j=0; j<n_bao; ++j ) {
                double delta_i = bao_data[i] - bao_prediction[i];
                double delta_j = bao_data[j] - bao_prediction[j];
                double icov_ij = bao_icov[i*n_bao+j];
                chi2 += delta_i*icov_ij*delta_j;
            }
        }
    }

    // DEBUG_TIME_END(_time_);

    return chi2;
}


double Likelihood_BAO(	imcmc_double&	params,
                        double&			lndet,
                        double&			chisq,
                        void*			model,
                        void*			data,
                        istate&			state ) {

    // time_t _time_;
    // DEBUG_TIME_START(_time_);

    lndet = chisq = 0.0;

    CosmoTheory	*cosmo 	= static_cast<CosmoTheory*>(model);
    Data_BAO	*bao 	= static_cast<Data_BAO*>(data);

    double chisq_temp;
    for( int i=0; i<bao->bao_data_size; ++i ) {

        double *bao_prediction = new double[bao->BAO_List[i]->size];

        Compute_BAO_Vals(	cosmo,
                            bao->BAO_List[i]->z,
                            bao_prediction,
                            bao->BAO_List[i]->size,
                            bao->BAO_List[i]->type,
                            bao->use_Hu_fitting ); // this place can be updated .... !!!!!!!!!!!!!!!!!

        chisq_temp = Compute_BAO_Chi2(	bao->BAO_List[i]->val,
                                        bao_prediction,
                                        bao->BAO_List[i]->err,
                                        bao->BAO_List[i]->icov,
                                        bao->BAO_List[i]->size );

//debug
//         bao->BAO_List[i]->print_type();
//         for(int j=0; j<bao->BAO_List[i]->size; ++j ){
//             cout << "data: " << bao->BAO_List[i]->val[j] 
//                  << "\ttheory: " << bao_prediction[j] 
//                  << ", ratio = " << bao->BAO_List[i]->val[j]/bao_prediction[j] 
//                  << endl;
//         }

        // cout <<"\n";

        chisq += chisq_temp;
        delete[] bao_prediction;
    }
    // exit(0);

    // DEBUG_TIME_END(_time_);

    return lndet - 0.5*chisq;
}
