/*!
 * Conjugate Gradient File
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cg.h"

#define LS_ITER_MAX 10
#define REC_MAX 10
#define DX_TOL 1.e-4
#define DELTAX_TOL 1.e1

void strong_wolfe( struct cg_workspace *g, double *x, double alpha0);
void  bracket_search( struct cg_workspace *g, double *x,
		      int sense, int rec_count);
double interpolate_cubic( double f0, double df0, double x0,
			  double f1, double df1, double x1);
double cg_beta_PRP( struct cg_workspace *g);
double estimate_alpha0( struct cg_workspace *g, double smag_old);
int stopping_condition( struct cg_workspace *g);


/*!
 * command to allocate the memory needed for the conjugate gradient method
 */
struct cg_workspace* allocate_cg_workspace( int n)
{
  struct cg_workspace *g =(struct cg_workspace *) 
    malloc( sizeof( struct cg_workspace));
  g->n = n;
  g->dx0 = (double *) malloc( n*sizeof(double) );
  g->s = (double *) malloc( n*sizeof(double) );
  g->dxp = (double *) malloc( n*sizeof(double) );

  return( g);
}

/*!
 * function that free's the memory needed for the conjugate gradient method
 */
void free_cg_workspace( struct cg_workspace *g)
{
  free( g->dxp);
  free( g->s);
  free( g->dx0);
  free( g);
}

/*!
 * Conjugate Gradient Algorithm
 * Taken from the wikipedia page
 */
void conjugate_gradient( struct cg_workspace *g, double *x)
{
  double (*h)( int n, const double *x, double *dx, void *p) = g->h;
  int n = g->n;
  void *p = g->params;

  int i, iter = 0, stopping_status;
  double alpha0, beta, dir, smag_temp;

  /* Evaluate the function at the current point */
  for( i=0; i<n; i++){
    (g->dx0)[i] = 0.;
  }
  g->h0 = h( n, x, g->dx0, p);
  for( i=0; i<n; i++){
    (g->dxp)[i] = (g->dx0)[i];
    (g->s)[i] = -(g->dx0)[i];
  }


  /* Get magnitude of the search direction vector */
  g->smag = 0.;
  for( i=0; i<n; i++){
    g->smag += (g->s)[i]*(g->s)[i];
  }
  g->smag = sqrt( g->smag);
  /* Make a first guess at alpha and do the line search */
  alpha0 = 1./(g->smag);
  strong_wolfe( g, x, alpha0);

  /* Evaluate the stopping condition */
  stopping_status = stopping_condition( g);

  while( !stopping_status ){
    /* evaluate the new search direction */
    beta = cg_beta_PRP( g);
    for( i=0; i<n; i++){
      (g->s)[i] = -(g->dxp)[i] + beta*(g->s)[i];
      (g->dx0)[i] = (g->dxp)[i];
    }
    /* check the reset condition for the CG method */
    dir = 0.;
    for( i=0; i<n; i++){
      dir += (g->s)[i]*(g->dx0)[i]; 
    }
    if( dir > 0. || iter > n ){
      iter = 0;
      for( i=0; i<n; i++){
	(g->s)[i] = -(g->dx0)[i];
      }
    }
    iter++;

    /* Get magnitude of the search direction vector */
    smag_temp = 0.;
    for( i=0; i<n; i++){
      smag_temp += (g->s)[i]*(g->s)[i];
    }
    smag_temp = sqrt( smag_temp);
    /* Estimate initial guess at alpha and do the line search */
    alpha0 = estimate_alpha0( g, smag_temp);
    g->smag = smag_temp;
    strong_wolfe( g, x, alpha0);
    
    /* Evaluate stopping condition */
    stopping_status = stopping_condition( g);
 }
}

/*!
 * Line search algorithm that satisfies the strong Wolfe conditions
 * taken from algorith 3.5 in .. and ..
 */
void strong_wolfe( struct cg_workspace *g, double *x, double alpha0)
{
  double (*h)( int n, const double *x, double *dx, void *p) = g->h;
  int n = g->n;
  void *p = g->params;

  double c1 = 1.e-4;
  double c2 = 0.1;
  int i, iter = 0;
  g->alpham = 0.;
  g->alphap = alpha0;

  /* Get f and df */
  g->hm = g->h0;
  g->dh0 = 0.;
  for( i=0; i<n; i++){
    g->dh0 += (g->dx0)[i]*(g->s)[i];
  }
  g->dh0 = g->dh0/(g->smag);
  g->dhm = g->dh0;

  do{
    for( i=0; i<n; i++){
      x[i] = x[i] + ((g->alphap)-(g->alpham))*(g->s)[i];
      (g->dxp)[i] = 0.;
    }
    g->hp = h( n, x, (g->dxp), p);
    g->dhp = 0.;
    for( i=0; i<n; i++){
      g->dhp += (g->dxp)[i]*(g->s)[i];
    }
    g->dhp = g->dhp/(g->smag);

    if( g->hp > (g->h0)+c1*(g->alphap)*(g->dh0) || 
	( g->hp >= g->hm && iter > 0) ){
      bracket_search( g, x, +1, 0);
      break;
    }
    if( fabs( g->dhp) <= -c2*(g->dh0) ){
      break;
    }
    if( g->dhp >= 0. ){
      bracket_search( g, x, -1, 0);
      break;
    }
    g->alpham = g->alphap;
    g->hm = g->hp;
    g->dhm = g->dhp;
    g->alphap = 2.*(g->alphap);
    iter++;
  }while( iter < LS_ITER_MAX );  
}

/*!
 * This is the bracketed search algorith
 * taken from algorith 3.6 in .. and ..
 * rewritten in recursive form
 */
void  bracket_search( struct cg_workspace *g, double *x,
		      int sense, int rec_count)
{
  double (*h)( int n, const double *x, double *dx, void *p) = g->h;
  int n = g->n;
  void *p = g->params;

  double c1 = 1.e-4;
  double c2 = 0.1;
  double fj, dfj, alphaj, flo;
  int i, success = 0, new_sense;

  fprintf( stdout, "range: %1.3g\n",
	   fabs( ((g->alpham)-(g->alphap))*(g->smag)) );

  alphaj = interpolate_cubic( g->hm, g->dhm, g->alpham,
  			      g->hp, g->dhp, g->alphap);
/*  alphaj = 0.5*((g->alphap)+(g->alpham)); */

  for( i=0; i<n; i++){
    x[i] = x[i] + (alphaj-(g->alphap))*(g->s)[i];
    (g->dxp)[i] = 0.;
  }
  fj = h( n, x, (g->dxp), p);
  dfj = 0.;
  for( i=0; i<n; i++){
    dfj += (g->dxp)[i]*(g->s)[i];
  }
  dfj = dfj/(g->smag);
  
  if( sense == +1 ){
    flo = g->hm;
  }else if( sense == -1){
    flo = g->hp;
  }

  if( fj > (g->h0) + c1*alphaj*(g->dh0) ||
      fj >= flo ){
    if( sense == -1 ){
      g->alpham = g->alphap;
      g->hm = g->hp;
      g->dhm = g->dhp;
    }
    (g->alphap) = alphaj;
    (g->hp) = fj;
    (g->dhp) = dfj;
    new_sense = +1;
  }else{
    if( fabs( dfj) <= -c2*(g->dh0) ){
      success = 1;
    }
    if( dfj*((g->alpham)-(g->alphap)) >= 0. ){
      g->alpham = g->alphap;
      g->hm = g->hp;
      g->dhm = g->dhp;
    }
    (g->alphap) = alphaj;
    (g->hp) = fj;
    (g->dhp) = dfj;
    new_sense = -1; 
  }
  if( success != 1 && rec_count < REC_MAX ){
    bracket_search( g, x, new_sense, rec_count+1);
  }
}

/*!
 * A truncated cubic interpolation.  This little function calculates
 * the location of the minimum of a cubic polynomial that has its
 * values and derivatives given at the end points. It the calculated
 * minimum lies outside of the center 3/4 of the interval then an
 * estimated position at either at 1/8th or 7/8th of the intevale is
 * given.
 */
double interpolate_cubic( double f0, double df0, double x0,
			  double f1, double df1, double x1)
{
  double d1, d2, sign, coef;
  double output;

  d1 = df0+df1-3.*(f0-f1)/(x0-x1);
  sign = ( (x1-x0) >= 0. ? 1. : -1.);
  d2 = sign*sqrt( d1*d1-df0*df1);
  coef = (df1+d2-d1)/(df1-df0+2.*d2);
  
  coef = ( coef > 0.125 ? coef : 0.125 );
  coef = ( coef < 0.875 ? coef : 0.875 );
  
  output = x1-(x1-x0)*coef;

  return( output );
}

double cg_beta_PRP( struct cg_workspace *g)
{
  int n = g->n;
  int i;
  double num, denom, beta;

  num = 0.;
  denom = 0.;
  for( i=0; i<n; i++){
    num += (g->dxp)[i]*((g->dxp)[i]-(g->dx0)[i]);
    denom += (g->dx0)[i]*(g->dx0)[i];
  }
  beta = num/denom;
  return ( beta > 0. ? beta : 0. );
}

double estimate_alpha0( struct cg_workspace *g, double smag_old)
{
  return( (g->alphap)*(g->smag)/smag_old);
}

int stopping_condition( struct cg_workspace *g)
{
  int n = g->n;
  int i;
  double deltax_linf_norm = fabs( (g->alphap)*(g->s)[0] );
  double dhdx_linf_norm = fabs( (g->dxp)[0]);

  for( i=1; i<n; i++){
    deltax_linf_norm = ( fabs( (g->alphap)*(g->s)[i]) > deltax_linf_norm ? 
			 fabs( (g->alphap)*(g->s)[i]) : deltax_linf_norm );
    dhdx_linf_norm = ( fabs( (g->dxp)[i]) > dhdx_linf_norm ? 
		       fabs( (g->dxp)[i]) : dhdx_linf_norm );
  }
  fprintf( stdout, " %1.3g %1.3g\n", deltax_linf_norm, dhdx_linf_norm);

  return (dhdx_linf_norm < DX_TOL && deltax_linf_norm < DELTAX_TOL);
}
