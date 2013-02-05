/*!
 *********************************************************************
 * Simple Non-linear Conjugate Gradient Optimization
 **********************************************************************
 *
 * This file implements my own version of a non-linear conjugate
 * gradient optimization routine. I wrote it using wikipedia as the
 * only reference. The method returns when either the change in the
 * vector is less than some tolerance, or the system has reaches a
 * maximum interation. Both the tollerance and iterations can be set
 * using the define macros.
 *
 * The beta in the conjugate gradient method is chosen as the max( 0,
 * beta_PR), where the beta_PR is fromthe Polak-Ribiere formula. 
 *
 * The line search is from a CS class in at U Wisconson, who had the
 * lecture notes on the web. The relevant lecture notes are included
 * in the something or other.
 */
#include <stdlib.h>
#include <math.h>

#include <stdio.h>

#define TOLERANCE 1.e-4
#define ITER_MAX 200
#define ALPHA_INIT 1.
#define ALPHA_MAX 1.e100

struct cg_workspace{
  double (*h)( int n, const double *x, double *dx, void *p);
  int n;
  void *params;

  double val;
  double alpha;

  double *dx;
  double *dx_old;
  double *s;
  double *s_old;

  double *x_alpha;
  double *dx_alpha;
};

double ls_function( struct cg_workspace *g, 
		    const double *x, const double *p 
		    double alpha, double *df)
{
  double f;
  int i;
  
  for( i=0; i<n; i++){
    g->x_alpha[i] = x[i] + alpha*p[i];
    g->dx_alpha[i] = 0.;
  }

  f = g->h( g->n, g->x_alpha, g->dx_alpha, g->params);
  *df = 0.;
  for( i=0; i<n; i++){
    *df += p[i]*dx_alpha[i];
  }

  return( f);
}

double line_search( struct cg_workspace *g,
		    const double *x, const double *p, double alpha_0)
{
  double f0, df0;
  double alpha_lo, f_lo, df_lo;
  double alpha_hi, f_hi, df_hi;

  alpha_lo = 0.;
  f0 = ls_function( g, x, p, 0., &df0);
  f_lo = f_0;
  df_lo = df_0;

  alpha_hi = alpha_0;
  f_hi = ls_function( g, x, p, alpha_hi, &df_hi);

}

struct cg_workspace* 
allocate_cg_workspace( double (*h)(int n, const double *x, 
				   double *dx, void *p),
		       int n, void *params)
{
  struct cg_workspace *g = (struct cg_workspace*) 
    malloc( sizeof(struct cg_workspace) );

  g->h = h;
  g->n = n;
  g->params = params;

  g->dx = (double *) malloc( n*sizeof(double));
  g->dx_old = (double *) malloc( n*sizeof(double));
  g->s = (double *) malloc( n*sizeof(double));
  g->s_old = (double *) malloc( n*sizeof(double));
  g->x_alpha = (double *) malloc( n*sizeof(double));
  g->dx_alpha = (double *) malloc( n*sizeof(double));

  return( g);
}

void free_cg_workspace( struct cg_workspace *g)
{
  free( g->dx);
  free( g->dx_old);
  free( g->s);
  free( g->s_old);
  free( g->x_alpha);
  free( g->dx_alpha);

  free( g);
}

double mag_array( int n, const double *x)
{
  int i;
  double mag = 0;
  
  for( i=0; i<n; i++){
    mag += x[i]*x[i];
  }
  
  return( sqrt(mag) );
}

double interpolate_alpha( double alpha_lo, double alpha_hi){
  return( (alpha_lo + alpha_hi)/2.);
}

double my_line_search( struct ls_workspace *s , double alpha_0){
  double alpha_lo = 0;
  double alpha_hi = alpha_0;
  



double zoom( struct cg_workspace *g, double *x,  
	     double phi_lo_i, double alpha_lo_i, 
	     double phi_hi_i, double alpha_hi_i)
{
  int i, iter=0;
  double phi_j, phip_j, phi_0, phip_0, phi_lo, phi_hi;
  double alpha, alpha_j, alpha_lo, alpha_hi;
  double c1 = 1.e-4, c2 = 0.9;

  phi_lo = phi_lo_i;
  alpha_lo = alpha_lo_i;
  phi_hi = phi_hi_i;
  alpha_hi = alpha_hi_i;

  phi_0 = g->val;
  phip_0 = 0;
  for( i=0; i<g->n; i++){
    phip_0 += -(g->s)[i]*(g->dx)[i];
  }
  
  while( iter < ITER_MAX){
    alpha_j = interpolate_alpha( alpha_lo, alpha_hi);

    for( i=0; i<g->n; i++){
      (g->x_alpha)[i] = x[i] + alpha_j*(g->s[i]);
      (g->dx_alpha)[i] = 0;
    }
    phi_j = g->h( g->n, g->x_alpha, g->dx_alpha, g->params);
    phip_j = 0;
    for( i=0; i<g->n; i++){
      phip_j += (g->s)[i]*(g->dx_alpha)[i];
    }

    if( ( phi_j > phi_0 + c1*alpha_j*phip_0) ||
	( phi_j >= phi_lo) ){
      alpha_hi = alpha_j;
      phi_hi = phi_j;
    }else{
      if( fabs(phip_j) <= -c2*phip_0 ){
	alpha = alpha_j;
	break;
      }
      if( phip_j*(alpha_hi-alpha_lo) >= 0. ){
	alpha_hi = alpha_lo;
	phi_hi = phi_lo;
      }
      alpha_lo = alpha_j;
      phi_lo = phi_j;
    }
  }

  return( alpha);
}

double line_search( struct cg_workspace *g, double *x, double alpha_1)
{
  int i, iter;
  double alpha, alpha_i, alpha_i_old;
  double phi_0, phi_i, phi_i_old;
  double phip_0, phip_i;
  double c1 = 1.e-4, c2 = 0.9;

  alpha_i_old = 0;
  alpha_i = alpha_1;

  phi_0 = g->val;
  phip_0 = 0;
  for( i=0; i<g->n; i++){
    phip_0 += -(g->s)[i]*(g->dx)[i];
  }
  phi_i = phi_0;

  iter = 0;
  while( iter < ITER_MAX ){

    phi_i_old = phi_i;
    for( i=0; i<g->n; i++){
      (g->x_alpha)[i] = x[i] + alpha_i*(g->s[i]);
      (g->dx_alpha)[i] = 0;
    }
    phi_i = g->h( g->n, g->x_alpha, g->dx_alpha, g->params);
    phip_i = 0;
    for( i=0; i<g->n; i++){
      phip_i += (g->s)[i]*(g->dx_alpha)[i];
    }
 
    if( ( phi_i > phi_0 + c1*alpha_i*phip_0) ||
	( phi_i >= phi_i_old && iter > 0) ){
      alpha = zoom( g, x, phi_i_old, alpha_i_old, 
		    phi_i, alpha_i);
      break;
    }else if( fabs(phip_i) <= -c2*phip_0 ){
      alpha = alpha_i;
      break;
    }else if( phip_i >= 0 ){
      alpha = zoom( g, x, phi_i, alpha_i, 
		    phi_i_old, alpha_i_old);
      break;
    }else{
      alpha_i_old = alpha_i;
      alpha_i = 2*alpha_i;
      iter++;
    }
  }

  return( alpha);
}

double get_beta( int n, const double *dx, 
		  const double *dx_old, const double *s)
{
  int i;
  double num, denom, beta;

  num=0;
  denom=0;
  for( i=0; i<n; i++){
    num += dx[i]*(dx[i]-dx_old[i]);
    denom += dx_old[i]*dx_old[i];
  }
  beta = num/denom;
  
  return( ( beta > 0 ) ? beta: 0);
}

double get_alpha_0( struct cg_workspace *g){
  int i;
  double num = 0., denom = 0.;

  for( i=0; i<g->n; i++){
    num += (g->s_old)[i]*(g->dx_old)[i];
    denom += (g->s)[i]*(g->dx)[i];
  }
 
  return( g->alpha*num/denom);
}

void cg_initialize( struct cg_workspace *g, double *x)
{
  int i;
  double alpha_0;
  
  for( i=0; i<g->n; i++){
    (g->dx)[i] = 0.;
  }
  g->val = g->h( g->n, x, g->dx, g->params);
  for( i=0; i<g->n; i++){
    (g->dx)[i] = -(g->dx)[i];
    (g->s)[i] = (g->dx)[i];
  }

  alpha_0 = ALPHA_INIT/mag_array( g->n, g->dx);
  g->alpha = line_search( g, x, alpha_0);
  for( i=0; i<g->n; i++){
    x[i] += g->alpha*(g->s)[i];
  }
}

void cg_step( struct cg_workspace *g, double *x)
{
  int i;
  double beta, alpha_0;

  for( i=0; i<g->n; i++){
    (g->dx_old)[i] = (g->dx)[i];
    (g->dx)[i] = 0.;
  }
  g->val = g->h( g->n, x, g->dx, g->params);
  for( i=0; i<g->n; i++){
    (g->dx)[i] = -(g->dx)[i];
  }

  beta = get_beta( g->n, g->dx, g->dx_old, g->s);
  for( i=0; i<g->n; i++){
    (g->s_old)[i] = (g->s)[i];
    (g->s)[i] = (g->dx)[i] + beta*(g->s)[i];
  }

  alpha_0 = get_alpha_0( g);
  g->alpha = line_search( g, x, alpha_0);

  for( i=0; i<g->n; i++){
    x[i] += g->alpha*(g->s)[i];
  }
}

void conjugate_gradient( struct cg_workspace *g, double *x)
{
  int iter=0, test1, test2;

  cg_initialize( g, x);
  test1 = (g->alpha*mag_array( g->n, g->s) > TOLERANCE);
  test2 = (iter< ITER_MAX);
  while( test1 && test2 ){
    cg_step( g, x);
    test1 = (g->alpha*mag_array( g->n, g->s) > TOLERANCE);
    test2 = (iter< ITER_MAX);
  }
}



