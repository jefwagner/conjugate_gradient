/*!*********************************************************************
 * \file cg_test.c
 * \author Jef Wagner (jefwagner@gmail.com)
 ***********************************************************************
 * This file contains the unit test for the functions used in the
 * conjugate gradient file.
 ***********************************************************************
 */

#include <stdio.h>
#include <cg.c>

/*!
 * Quadratic test function parameter struct.
 */
typedef struct{ double *a, *c, fmin;} quad_params;

/*!
 * Simple quadratic test function.
 *
 * This simple test function for optimization, 

 * \f[ f({\bf x\}) = \sum_i a_i (x_i-c_i)^2 + f_{\text{min}}, \f], 

 * for all \f$a\f$ positive, has minimum \f$ f({\bf x}) =
 * f_{\text{min}}\f$ at \f${\bf x} = {\bf c}\f$.
 */
double quad( unsigned int n, const double *x, 
	     double *dfdx, void *p){
  int i;
  double f;
  quad_params *qp = (quad_params*) p;
  double *a = qp->a;
  double *c = qp->c;
  f = qp->fmin;
  for( i=0; i<n; i++){
    f += a[i]*(x[i]-c[i])*(x[i]-c[i]);
    dfdx[i] = 2.*a[i]*(x[i]-c[i]);
  }
  return f;
}

/*!
 * opt_fn_eval test.
 *
 * This function test to see if the opt_fn_eval function works as
 * advertized.
 *
 * We use the simple quadratic test function with \f$ n = 10 \f$, \f$
 * {\bf a} = \{ 1, 1, \ldots, 1 \} \f$, and \f$ {\bf c} = 0 \f$. In
 * this configuration the function gives the square of the Euclidean
 * size of the argument \f$ f({\bf x}) = |{\bf x}|^2 \f$, and the
 * gradient is equal to the argument \f$
 * \frac{\mathrm{d}f}{\mathrm{d}{\bf x}} = {\bf x} \f$. In order to
 * test this, we succesively test

 * - \f$ {\bf x} = \{ 1, 0, 0, 0, \ldots, 0 \},
 * - \f$ {\bf x} = \{ 1, 1, 0, 0, \ldots, 0 \},
 * - \f$ {\bf x} = \{ 1, 1, 1, 0, \ldots, 0 \},

 * until all elements are 1. We test to see if the value of the
 * function gives the number of 1 elements, and see if the gradient
 * matches the argument.
 */
void opt_fn_eval_test(){
  int i, j, n=10;
  int status;
  double f, *x, *dfdx;
  quad_params qp;
  opt_fn of;
  
  fprintf( stdout, " opt_func_eval: ");

  /* Allocate memory for arrays */
  x = (double *) malloc( n*sizeof( double));
  dfdx = (double *) malloc( n*sizeof( double));
  qp.a = (double *) malloc( n*sizeof( double));
  qp.c = (double *) malloc( n*sizeof( double));

  /* Fill parameter arrays */
  for( i=0; i<n; i++){
    qp.a[i] = 1.;
    qp.c[i] = 0.;
  }
  qp.fmin = 0.;

  /* Assign members to opt_fn structure */
  of.n = n;
  of.f = (objective_fn) &quad;
  of.p = (void*) &qp;

  /* Initialize data and status */
  for( i=0; i<n; i++){
    x[i] = 0.;
  }
  status = 0;
  /* Loop through all dimensions */
  for( i=0; i<n; i++){
    /* Change ith component to 1 */
    x[i] = 1.;
    f = opt_fn_eval( x, dfdx, &of);
    /* test that f(x) = i+1 */
    if( fabs( f-(i+1)) > 1.e-6 ){
      status++;
    }
    /* test that dfdx = 2x */
    for( j=0; j<n; j++){
      if( fabs( 2.*x[i]-dfdx[i]) > 1.e-6 ){
	status++;
      }
    }
  }
  /* pass if status = 0, otherwise fail */
  if( status == 0 ){
    fprintf( stdout, "pass\n" );
  }else{
    fprintf( stdout, "fail\n" );
  }

  /* free the dynamically allocated memory */
  free( qp.c);
  free( qp.a);
  free( dfdx);
  free( x);
}

void lin_fn_eval_test(){
  double a[2], c[2], fmin;
  quad_params qp;
  opt_fn of;
  double x[2], dfdx[2], s[2];
  lin_fn lf;
  
  a[0] = ;
  a[1] = ;
  c[0] = ;
  c[1] = ;
  qp.a = a;
  qp.c = c;
  qp.fmin = ;
  
  of.n = 2;
  of.f = (objective_func) &quad;
  of.p = (void*) &qp;

  x[0] = ;
  x[1] = ;
  dfdx[0] = ;
  dfdx[1] = ;

  lf.of = of;
  lf.x = x;
  lf.dfdx = dfdx;
  lf.s = s;

}

/*
  for( i=0; i<n; i++){
    qp.a[i] = 0.5 + (1.*rand())/RAND_MAX;
    qp.c[i] = (2.*rand())/RAND_MAX -1.;
  }
  qp.fmin = (2.*rand())/RAND_MAX -1.;
*/

int main(){
  opt_fn_eval_test();
  lin_fn_eval_test();
  return 0;
}
