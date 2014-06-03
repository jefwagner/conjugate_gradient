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

/*!
 * 1-D quadratic test function.
 */
double quad_1d( unsigned int n, const double *x, 
		double *dfdx, void *p){
  dfdx[0] = 2.*x[0];
  return x[0]*x[0];
}
/*!
 * 1-D gaussian test function
 */
double gaussian_1d( unsigned int n, const double *x, 
		    double *dfdx, void *p){
  double f = exp( -0.5*x[0]*x[0]);
  dfdx[0] = x[0]*f;
  return -f;
}
/*!
 * 1-D quartic double well test function
 */
double quartic_1d( unsigned int n, const double *x,
		   double *dfdx, void *p){
  double x2 = x[0]*x[0];
  dfdx[0] = 4.*x[0]*x2-2.*x[0];
  return x2*x2-x2;
}

/*!
 * rand_double utility function
 */
double rand_double( double low, double hi){
  return (hi-low)*rand()/RAND_MAX + low;
}

/*!
 * sw_bracket_search test
 *
 * This function uses the three simple 1-D test functions to test the
 * bracket search algorithm. The starting points for the bracketed
 * search chosen randomly. For the quadratic we also randomly choose
 * between the order of the arguments. The test succeeds if the final
 * points obey the strong wolfe conditions.
 */
void sw_bracket_search_test(){
  int status, i;
  double c1=C1, c2=C2;
  double x, f, dfdx;
  double x0, ap, am;
  double f0, df0, fm, dfm, fp, dfp;
  opt_fn of;
  lin_fn lf;
  double s, smag;

  fprintf( stdout, " sw_bracket_search: ");
  status = 0;

  /* Test the quadratic function*/
  of.n = 1;
  of.f = (objective_fn) &quad_1d;
  lf.of = of;
  
  for( i=0; i<10; i++){
    /* Choose the initial point randomly */
    x0 = rand_double( 1., 6.);

    /* Evaluate the function at x0 */
    f0 = opt_fn_eval( &x0, &df0, &of);

    /* Set the info for the linear helper struct */
    x = x0;
    s = -df0;
    smag = fabs( df0);
    lf.x = &x;
    lf.dfdx = &dfdx;
    lf.s = &s;
    lf.smag = smag;
    lf.a_prev = 0.;
    
    /* Set initial data to check wolf conditions*/
    f0 = lin_fn_eval( 0, &df0, &lf);
    /* Always use one bracket point as our initial point */
    am = 0;
    fm = f0;
    dfm = df0;
    /* with 50% chance */
    if( rand() % 2 == 0 ){
      /* Choose a new point that is "larger" */
      ap = rand_double( 2*(1-c1)*x0, 4*(1-c1)*x0);
      fp = lin_fn_eval( ap, &dfp, &lf);
      
      f = sw_bracket_search( f0, df0, am, fm, dfm, ap, fp, dfp, &lf);
    }else{
      /* Choose a new point that is "smaller" */
      ap = rand_double( (1+c2)*x0, 2*(1-c1)*x0);
      fp = lin_fn_eval( ap, &dfp, &lf);
      
      f = sw_bracket_search( f0, df0, ap, fp, dfp, am, fm, dfm, &lf);
    }
    /* Check if the result satisfies the wolfe conditions */
    if( f > f0 + c1*(x0-x)*df0 || fabs( dfdx) > -c2*df0 ){
      status += 1;
    } 
  }

  /* Test the gaussian function*/
  of.n = 1;
  of.f = (objective_fn) &gaussian_1d;
  lf.of = of;
  
  for( i=0; i<10; i++){
    /* Choose the initial point randomly */
    x0 = rand_double( 1., 3.);

    /* Evaluate the function at x0 */
    f0 = opt_fn_eval( &x0, &df0, &of);

    /* Set the info for the linear helper struct */
    x = x0;
    s = -df0;
    smag = fabs( df0);
    lf.x = &x;
    lf.dfdx = &dfdx;
    lf.s = &s;
    lf.smag = smag;
    lf.a_prev = 0.;
    
    /* Set initial data to check wolf conditions*/
    f0 = lin_fn_eval( 0, &df0, &lf);
    /* Always use one bracket point as our initial point */
    am = 0;
    fm = f0;
    dfm = df0;
    /* Choose a new point that is "larger" */
    ap = rand_double( 2*(1-c1)*x0 , 4*(1-c1)*x0);
    fp = lin_fn_eval( ap, &dfp, &lf);
      
    f = sw_bracket_search( f0, df0, am, fm, dfm, ap, fp, dfp, &lf);
    /* Check if the result satisfies the wolfe conditions */
    if( f > f0 + c1*(x0-x)*df0 || fabs( dfdx) > -c2*df0 ){
      status += 1;
    } 
  }

  /* Test the quartic function*/
  of.n = 1;
  of.f = (objective_fn) &quartic_1d;
  lf.of = of;
  
  for( i=0; i<10; i++){
    /* Choose the initial point randomly */
    x0 = rand_double( 1., 3.);

    /* Evaluate the function at x0 */
    f0 = opt_fn_eval( &x0, &df0, &of);

    /* Set the info for the linear helper struct */
    x = x0;
    s = -df0;
    smag = fabs( df0);
    lf.x = &x;
    lf.dfdx = &dfdx;
    lf.s = &s;
    lf.smag = smag;
    lf.a_prev = 0.;
    
    /* Set initial data to check wolf conditions*/
    f0 = lin_fn_eval( 0, &df0, &lf);
    /* Always use one bracket point as our initial point */
    am = 0;
    fm = f0;
    dfm = df0;
    /* Choose a new point that is "larger" */
    ap = rand_double( 2*(1-c1)*x0 , 4*(1-c1)*x0);
    fp = lin_fn_eval( ap, &dfp, &lf);
      
    f = sw_bracket_search( f0, df0, am, fm, dfm, ap, fp, dfp, &lf);
    /* Check if the result satisfies the wolfe conditions */
    if( f > f0 + c1*(x0-x)*df0 || fabs( dfdx) > -c2*df0 ){
      status += 1;
    } 
  }
  
  /* pass if status = 0, otherwise fail */
  if( status == 0 ){
    fprintf( stdout, "pass\n" );
  }else{
    fprintf( stdout, "fail\n" );
  }
}

/*!
 * sw_line_search test.
 *
 * This function test the strong wolfe line search proceedure using
 * the same 1-d test functions as the bracketed search function. The
 * starting points are randomly choosen. The test pass if the final
 * points obey the strong wolfe conditions.
 */
void sw_line_search_test(){
  int status, i;
  double c1=C1, c2=C2;
  double x, f, dfdx;
  double x0, f0, df0;
  opt_fn of;
  lin_fn lf;
  double s, smag;

  fprintf( stdout, " sw_line_search: ");
  status = 0;

  /* Test the quadratic function*/
  of.n = 1;
  of.f = (objective_fn) &quad_1d;
  lf.of = of;
  
  for( i=0; i<10; i++){
    /* Choose the initial point randomly */
    x0 = rand_double( 1., 6.);

    /* Evaluate the function at x0 */
    f0 = opt_fn_eval( &x0, &df0, &of);
    /* Change the sign of derivative, because the search direction is
       negative */
    df0 = -df0;

    /* Set the info for the linear helper struct */
    x = x0;
    s = df0;
    smag = fabs( df0);
    lf.x = &x;
    lf.dfdx = &dfdx;
    lf.s = &s;
    lf.smag = smag;
    lf.a_prev = 0.;
    
    f = sw_line_search( f0, &lf);
    if( f > f0 + c1*(x0-x)*df0 || fabs( dfdx) > -c2*df0 ){
      status += 1;
    } 
  }

  /* Test the gaussian function*/
  of.n = 1;
  of.f = (objective_fn) &gaussian_1d;
  lf.of = of;
  
  for( i=0; i<10; i++){
    /* Choose the initial point randomly */
    x0 = rand_double( 1., 3.);

    /* Evaluate the function at x0 */
    f0 = opt_fn_eval( &x0, &df0, &of);
    /* Change the sign of derivative, because the search direction is
       negative */
    df0 = -df0;

    /* Set the info for the linear helper struct */
    x = x0;
    s = df0;
    smag = fabs( df0);
    lf.x = &x;
    lf.dfdx = &dfdx;
    lf.s = &s;
    lf.smag = smag;
    lf.a_prev = 0.;
    
    f = sw_line_search( f0, &lf);
    /* Check if the result satisfies the wolfe conditions */
    if( f > f0 + c1*(x0-x)*df0 || fabs( dfdx) > -c2*df0 ){
      status += 1;
    } 
  }

  /* Test the quartic function*/
  of.n = 1;
  of.f = (objective_fn) &quartic_1d;
  lf.of = of;
  
  for( i=0; i<10; i++){
    /* Choose the initial point randomly */
    x0 = rand_double( 1., 3.);

    /* Evaluate the function at x0 */
    f0 = opt_fn_eval( &x0, &df0, &of);
    /* Change the sign of derivative, because the search direction is
       negative */
    df0 = -df0;

    /* Set the info for the linear helper struct */
    x = x0;
    s = df0;
    smag = fabs( df0);
    lf.x = &x;
    lf.dfdx = &dfdx;
    lf.s = &s;
    lf.smag = smag;
    lf.a_prev = 0.;
    
    f = sw_line_search( f0, &lf);
    /* Check if the result satisfies the wolfe conditions */
    if( f > f0 + c1*(x0-x)*df0 || fabs( dfdx) > -c2*df0 ){
      status += 1;
    } 
  }
  
  /* pass if status = 0, otherwise fail */
  if( status == 0 ){
    fprintf( stdout, "pass\n" );
  }else{
    fprintf( stdout, "fail\n" );
  }
}

/*!
 * nlcg_optimize test.
 *
 * This test three functions: nlcg_malloc, nlcg_set, and
 * nlcg_optimize.  It allocates a workspace for a function of up to
 * 100 variables, and confirms that it has a valid pointer. It then
 * sets the size of the system from 10 up to 110 variables in steps of
 * 10. For 10 through 100 we confirm that no error is returned, for
 * 110 we confirm that an error is returned. We randomly populate a
 * multidimensional quadratic function, and optimize. We confirm that
 * the value, position, and derivatives are all within tolerance.
 */
void nlcg_optimize_test(){
  int i, j, status, err;
  quad_params qp;
  double x[110], a[110], c[110];
  nlcg_ws g; 
  double f, dx, dfdx;
  double f_tol = 1.e-4;
  double dx_tol = 1.e-4;
  double dfdx_tol = 1.e-8;
  
  status = 0;

  qp.a = a;
  qp.c = c;

  /* allocate up to a 100 length arrays */
  g = nlcg_malloc( 100);
  /* check if it allocated correctly */
  fprintf( stdout, " nlcg_malloc: ");
  if( g == NULL ){ 
    fprintf( stdout, "fail\n" );
    return;
  }
  fprintf( stdout, "pass\n" );

  /* for different lengths */
  for( i=10; i<=110; i+=10){
    /* set the length */
    err = nlcg_set( (objective_fn) &quad, i, &qp, g);
    /* make sure it worked for i <= 100 */
    if( err != NLCG_SUCCESS ){
      if( i==110){
	fprintf( stdout, " nlcg_set: pass\n" );
      }
      break;
    }
    /* make sure that didn't work for i=110 */
    if( i==110 ){
      fprintf( stdout, " nlcg_set: FAIL\n" );
      break;
    }
    /* randomly fill up the quadratic parameters */
    for( j=0; j<i; j++){
      a[j] = rand_double( 0.1, 2.1);
      c[j] = rand_double( -3., 3.);
      x[j] = rand_double( -3., 3.);
    }
    qp.fmin = rand_double( 0., 1.);
    /* optimize the quadratic function */
    f = nlcg_optimize( x, g);
    /* Check if the solution value is near the minimum value */
    if( fabs(f-qp.fmin) > f_tol ){
      status += 1;
    }
    /* Check if the solution is near the minimum */
    dx = 0.;
    for( j=0; j<i; j++){
      dx += (x[j]-c[j])*(x[j]-c[j]);
    }
    dx = sqrt( dx);
    if( dx > dx_tol ){
      status += 1;
    }
    /* Check if the derivative is near zero */
    dfdx = 0.;
    for( j=0; j<i; j++){
      dfdx = (g->lf.dfdx[j])*(g->lf.dfdx[j]);
    }
    dfdx = sqrt( dfdx);
    if( dfdx > dfdx_tol ){
      status += 1;
    }
  }

  fprintf( stdout, " nlcg_optimze: ");
  /* pass if status = 0, otherwise fail */
  if( status == 0 ){
    fprintf( stdout, "pass\n" );
  }else{
    fprintf( stdout, "fail\n" );
  }

  nlcg_free( g); 
}

int main(){
  opt_fn_eval_test();
  lin_fn_eval_test();
  sw_bracket_search_test();
  sw_line_search_test();
  nlcg_optimize_test();

  return 0;
}


