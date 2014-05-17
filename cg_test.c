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

/*!
 * lin_fn_eval test.
 *
 * Here we test to see if the linear function evaluation function
 * works as advertized. By definition, the linear function is

 * \f[ \phi(\alpha) = f(x_0+\alpha \hat{s}_x, y_0+\alpha \hat{s}_y), \]

 * where (x_0,y_0) is the starting point, and \f$\hat{{\bf s}}\f$ is a
 * normalized search direction. We use a simple two dimentional
 * quadratic function

 * \f[ f(x,y) = A x^2 + B y^2, \f]

 * Where A and B are choosen to be different numbers. We choose a
 * search direction in the (1,1) direction. where we know that the
 * linear function should give

 * \f[ \phi(\alpha) = \frac{A+B}{2}\alpha^2 + 
 *    \sqrt{2)(A x_0 + B y_0)\alpha + \big( A x_0^2 + B y_0^2\big). \f]

 * We will test several values of A and B, several magnitudes of the
 * search direction, and several starting points.
 */
void lin_fn_eval_test(){
  int i, j, status;
  double a[2], c[2];
  quad_params qp;
  opt_fn of;
  double x[2], dfdx[2], s[2], smag;
  lin_fn lf;
  double A, B, C;
  double alpha, f, tf, dfda, tdfda;
  
  fprintf( stdout, " lin_func_eval: ");

  status = 0;
  for( i=0; i<10; i++){

    a[0] = 5.*rand()/RAND_MAX;
    a[1] = 5.*rand()/RAND_MAX;
    c[0] = 0.;
    c[1] = 0.;
    qp.a = a;
    qp.c = c;
    qp.fmin = 0.;
  
    of.n = 2;
    of.f = (objective_fn) &quad;
    of.p = (void*) &qp;

    x[0] = 10.*rand()/RAND_MAX-5.;
    x[1] = 10.*rand()/RAND_MAX-5.;
    smag = exp( 6.*rand()/RAND_MAX-3.);
    s[0] = sqrt(0.5)*smag;
    s[1] = sqrt(0.5)*smag;

    lf.of = of;
    lf.x = x;
    lf.dfdx = dfdx;
    lf.s = s;
    lf.smag = smag;
    lf.a_prev = 0.;

    A = 0.5*(a[0]+a[1]);
    B = sqrt(2.)*(a[0]*x[0]+a[1]*x[1]);
    C = a[0]*x[0]*x[0]+a[1]*x[1]*x[1];
    
    for( j=0; j<10; j++){
      alpha = 5.*rand()/RAND_MAX;
      tf = A*alpha*alpha+B*alpha+C;
      tdfda = 2.*A*alpha+B;
      f = lin_fn_eval( alpha, &dfda, &lf);
      if( fabs( f-tf ) > 1.e-6 ){
	status += 1;
      }
      if( fabs( dfda-tdfda) > 1.e-6 ){
	status += 1;
      }
    }
  }
  if( status == 0 ){
    fprintf( stdout, "pass\n" );
  }else{
    fprintf( stdout, "fail\n" );
  }
}

double quad_1d( unsigned int n, const double *x, 
		double *dfdx, void *p){
  dfdx[0] = 2.*x[0];
  return x[0]*x[0];
}
double gaussian_1d( unsigned int n, const double *x, 
		    double *dfdx, void *p){
  double f = exp( -0.5*x[0]*x[0]);
  dfdx[0] = x[0]*f;
  return -f;
}
double quartic_1d( unsigned int n, const double *x,
		   double *dfdx, void *p){
  double x2 = x[0]*x[0];
  dfdx[0] = 3.*x[0]*x2-2.*x[0];
  return x2*x2-x2;
}

void sw_bracket_search_test(){
  
}



int main(){
  opt_fn_eval_test();
  lin_fn_eval_test();
  return 0;
}
