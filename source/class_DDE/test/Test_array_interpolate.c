#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _FAILURE_ -1
#define _SUCCESS_ 1
#define _ERRORMSGSIZE_ 1024
typedef char ErrorMsg[_ERRORMSGSIZE_];

int array_interpolate(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   double x,
		   int * last_index,
		   double * result,
		   int result_size, /** from 1 to n_columns */
		   ErrorMsg errmsg) {

  int inf,sup,mid,i;
  double weight;

  inf=0;
  sup=n_lines-1;

  if (*(array+inf*n_columns+index_x) < *(array+sup*n_columns+index_x)){
    if (x < *(array+inf*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array+inf*n_columns+index_x));
      return _FAILURE_;
    }
    if (x > *(array+sup*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array+sup*n_columns+index_x));
      return _FAILURE_;
    }
    while (sup-inf > 1) {
      mid=(int)(0.5*(inf+sup));
      if (x < *(array+mid*n_columns+index_x)) {sup=mid;}
      else {inf=mid;}
    }
  }
  else {
    if (x < *(array+sup*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array+sup*n_columns+index_x));
      return _FAILURE_;
    }
    if (x > *(array+inf*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array+inf*n_columns+index_x));
      return _FAILURE_;
    }
    while (sup-inf > 1) {
      mid=(int)(0.5*(inf+sup));
      if (x > *(array+mid*n_columns+index_x)) {sup=mid;}
      else {inf=mid;}
    }
  }

  *last_index = inf;

  weight=(x-*(array+inf*n_columns+index_x))/(*(array+sup*n_columns+index_x)-*(array+inf*n_columns+index_x));

/*    printf("x = %10.5f\t", x);*/

  for (i=0; i<result_size; i++){
    *(result+i) = *(array+inf*n_columns+i) * (1.-weight)
      + weight * *(array+sup*n_columns+i);
   
/*    printf("*(result+%d) = %10.5f\t",i,*(result+i));   */
  }
  
/*  printf("\n");*/

  *(result+index_x) = x;

  return _SUCCESS_;
}

int main(){

    int nLines = 21;
    int nColns = 2;
    double a =-3.14;
    double b = 3.14;
    double xi,yi;
    double *array = (double *)malloc(sizeof(double)*nLines*nColns);
    double *array_result = (double *)malloc(sizeof(double)*2);

    for( int i=0; i<nLines; i++ ){
        xi = a + i*(b-a)/(nLines-1.);
        yi = sin(xi);
        *(array + i*nColns + 0) = xi;
        *(array + i*nColns + 1) = yi;
        
/*        printf("xi = %10.5f\tyi = %10.5f\n",xi,yi);*/
    }
/*    exit(0);*/
    
    xi = a;
    int last_index;
    for( int i=0; i<100; i++ ){
        xi = a+i*(b-a)/(100.-1.0);
        array_interpolate( array, nColns, nLines, 0, xi, &last_index, array_result, 2, "");
        printf("%10.5f\t%10.5f\n",*(array_result+0),*(array_result+1));
    };

    
    free(array);
    free(array_result);
    
    return 0;
}
