#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "class.h"

#define index_DDE_z 0
#define index_DDE_w 1
#define index_DDE_w_over_1_and_z 2
#define index_DDE_int_w_over_1_and_z 3
#define index_DDE_splined_w 4
#define index_DDE_splined_w_over_1_and_z 5
#define index_DDE_splined_int_w_over_1_and_z 6
#define index_DDE_dervie_w_over_1_and_z 7

double w0 = -1.2;
double wa = 0.8;
double zmax = 1.5;

ErrorMsg errmsg;

int make_EOS_table( double **ET, double w0, double wa, double zmax, int n_lines, int n_columns ){

    double *EoS_table = (double *)malloc(sizeof(double)*n_lines*n_columns);

    for( int i=0; i<n_lines; i++ ){
        double z = i*zmax/(n_lines-1.0);
        double w = w0 + wa*z/(1.+z);
        *(EoS_table + i*n_columns + 0) = w;
        *(EoS_table + i*n_columns + 1) = w/(1.+z);
        *(EoS_table + i*n_columns + 2) = 0.;
        // *(EoS_table + i*n_columns + index_DDE_splined_w) = 0.;
        // *(EoS_table + i*n_columns + index_DDE_splined_w_over_1_and_z) = 0.;
        // *(EoS_table + i*n_columns + index_DDE_splined_int_w_over_1_and_z) = 0.;
    }

    *ET = EoS_table;

    return 0;
}

int make_EOS_spline_table( double **ETS, int n_lines, int n_columns ){

    double *EoS_spline_table = (double *)malloc(sizeof(double)*n_lines*n_columns);
    *ETS = EoS_spline_table;

    return 0;
}

int free_EOS_table( double *EoS_table ){
    if( EoS_table != NULL ){
        free(EoS_table);
        EoS_table = NULL;
    }
    return 0;
}

int free_EOS_spline_table( double *EoS_spline_table ){
    if( EoS_spline_table != NULL ){
        free(EoS_spline_table);
        EoS_spline_table = NULL;
    }
    return 0;
}

int main(){

    int n_lines = 20;
    int n_columns = 6;

    double *temp = (double *)malloc(sizeof(double)*n_lines*6);
    double *z_array = (double *)malloc(sizeof(double)*n_lines);
    
    double *EoS_table, *EoS_spline_table;
    make_EOS_table(&EoS_table,w0,wa,zmax,n_lines,3);
    make_EOS_spline_table(&EoS_spline_table,n_lines,3);

    for( int i=0; i<n_lines; i++){
        
        z_array[i] = i*zmax/(n_lines-1.);

        temp[i*n_columns+0] = EoS_table[i*3+0];
        temp[i*n_columns+1] = EoS_table[i*3+1];
        temp[i*n_columns+2] = EoS_table[i*3+2];

        // printf("EoS_table = %10.5f %10.5f %10.5f\n",
        //                     EoS_table[i*3+0],
        //                     EoS_table[i*3+1],
        //                     EoS_table[i*3+2]);
    }
    // exit(0);


/* step 1) initialize spline_array for w, w/1+z */
    class_call(array_spline_table_line_to_line(
                            z_array,
                            n_lines,
                            temp,
                            n_columns,
                            0,
                            3,
                            _SPLINE_EST_DERIV_,
                            errmsg),
                errmsg,
                errmsg);

    class_call(array_spline_table_line_to_line(
                            z_array,
                            n_lines,
                            temp,
                            n_columns,
                            1,
                            4,
                            _SPLINE_EST_DERIV_,
                            errmsg),
                errmsg,
                errmsg);

    class_call(array_integrate_spline_table_line_to_line(
                            z_array,
                            n_lines,
                            temp,
                            n_columns,
                            1,
                            4,
                            2,
                            errmsg),
                errmsg,
                errmsg);


    for( int i=0; i<n_lines; i++ ){
        if(i==0){
            EoS_table[i*3+2] = temp[i*n_columns+0];    
        }
        else
            EoS_table[i*3+2] = temp[i*n_columns+2]/log(1.+z_array[i]);
        // printf("z = %6.4f EoS_table[%2d*3+2] = %6.4f\n", z_array[i], i, EoS_table[i*3+2]);
    }

    class_call(array_spline_table_lines(
                            z_array,
                            n_lines,
                            EoS_table,
                            3,
                            EoS_spline_table,
                            _SPLINE_EST_DERIV_,
                            errmsg),
                errmsg,
                errmsg);

    // save interpolating results into a txt file
    FILE *fp = fopen("eos.txt","w");

    double z=0.0;
    double dz = 0.01;
    double w, w_over_1_and_z, int_w_over_1_and_z, splined_int_w_over_1_and_z, derive_w_over_1_and_z;
    double result[3] =  {0,0,0};
    int last_index;

    // while( z <= zmax ){
    //     class_call(array_interpolate_spline(
    //                         z_array,
    //                         n_lines,
    //                         EoS_table,
    //                         EoS_spline_table,
    //                         3,
    //                         z,
    //                         &last_index,
    //                         result,
    //                         3,
    //                         errmsg),
    //                 errmsg,
    //                 errmsg);

    //     fprintf(stdout,"%6.4f %6.4f %6.4f\t%6.4f\n",z,result[0],result[1],result[2]);
    //     z += dz;
    // };

    z=10;


    // int test = array_interpolate_spline(
    //                             z_array,
    //                             n_lines,
    //                             EoS_table,
    //                             EoS_spline_table,
    //                             3,
    //                             z,
    //                             &last_index,
    //                             result,
    //                             3,
    //                             errmsg);

    // if( test == _FAILURE_ ){
    //     fprintf(stdout,"%s\n",errmsg);
    //     exit(0);
    // }

    // class_test( z < 0.0 || z > zmax, ....)

    class_call( array_interpolate_spline(
                                z_array,
                                n_lines,
                                EoS_table,
                                EoS_spline_table,
                                3,
                                z,
                                &last_index,
                                result,
                                3,
                                errmsg),
                errmsg,
                errmsg);



    fprintf(stdout,"%6.4f %6.4f %6.4f\t%6.4f\n",z,result[0],result[1],result[2]);

    fclose(fp);

    free_EOS_table(EoS_table);
    free(EoS_spline_table);
    free(z_array);
    
    return 0;
}
