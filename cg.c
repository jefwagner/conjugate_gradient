/*!*********************************************************************
 * \file cg.c
 * \author Jef Wagner (jefwagner@gmail.com)
 ***********************************************************************
 * A very simple implementation of the non-linear conjugate gradient
 * method for the unconstraned optimization. The algorithm for the
 * non-linear conjugate gradient part was taken from the [wikipedia
 * article][NLCG]. The line search algorithms are algorithms 3.5 and
 * 3.6 taken from chapter 3 of [Numerical Optimization][LO] by Nocedal
 * and Write. The line search algorithm searches until the [strong
 * Wolfe conditions][SWC] is satisfied.
 ***********************************************************************
 * [NLCG] https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
 * [LO] https://www.amazon.com//dp/0387303030
 * [SWC] https://en.wikipedia.org/wiki/Wolfe_conditions
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "cg.h"

#define NLCG_ITER_MAX 100000
#define LS_ITER_MAX 20
#define BRACKET_ITER_MAX 40
#define ALPHA_0 1.
#define DFDX_TOL 1.e-8
#define C1 1.e-4
#define C2 0.1

/*
 * Optimization function helper structure.
 *
 * This structure holds the parameters that are needed to call the
 * objective function so that the function call can be
 * simplified. This is used in the `opt_fn_eval` function.
 */
typedef struct opt_fn_struct{
  objective_fn f; /*!< objective function */
  unsigned int n; /*!< size \f${\bf x}\f$ vector */
  void *p; /*!< additional parameters */
} opt_fn;

/*
 * Linear function helper structure
 *
 * This structure holds the parameters neccessary to evaluate a function

 * \f[ \phi(\alpha) = f({\bf x} + \alpha {\bf s}). \f]

 * This is used in the `lin_fn_eval` function, the `sw_line_search`
 * and `sw_bracket_search` functions.
 */
typedef struct lin_fn_struct{
  opt_fn of; /*!< optimization function helper struct */
  double *x; /*!< the search point vector */
  double *dfdx; /*!< the gradient vector */

  double *s; /*!< the search direction */
  double smag; /*!< the magnitude of the search direction */
  double a_prev; /*!< the last evaluated a value */
} lin_fn;

/*!
 * Non-linear conjugate gradient workspace.
 *
 * This function holds the neccessary variables and dynamically
 * allocated memory needed for the non-linear conjugate gradient
 * optimization. It is instantiated with `nlcg_malloc`, is freed with
 * `nlcg_free`, and it is used with `nlcg_set` and `nlcg_optimize`.
 */
struct nlcg_ws_struct{
  lin_fn lf; /*<! line search helper struct */
  double *dfdx_old; /*<! old gradient vector */
  unsigned int max_size; /*<! maximum size */
};


static double opt_fn_eval( const double *x, double *dfdx, const opt_fn *of);
static double lin_fn_eval( double a, double *dfda, lin_fn *lf);
static double sw_line_search( double f, lin_fn *lf);
static double sw_bracket_search( double f_0, double df_0,
				 double a_lo, double f_lo, double df_lo,
				 double a_hi, double f_hi, double df_hi, 
				 lin_fn *lf);
/*!
 * Non-linear conjugate gradient constructor.
 *
 * This function allocates the memory for the non-linear conjugate
 * gradient workspace. It returns a opaque pointer to the workspace.
 */
nlcg_ws nlcg_malloc( unsigned int max_size){
  nlcg_ws g = (nlcg_ws) malloc( sizeof(struct nlcg_ws_struct));
  if( g == NULL){ return NULL;}

  g->max_size = max_size;
  g->dfdx_old = (double*) malloc( max_size*sizeof(double));
  if( g->dfdx_old == NULL ){ return NULL;}
  g->lf.dfdx = (double*) malloc( max_size*sizeof(double));
  if( g->lf.dfdx == NULL ){ return NULL;}
  g->lf.s = (double*) malloc( max_size*sizeof(double));
  if( g->lf.dfdx == NULL ){ return NULL;}

  return g;
} 

/*!
 * Non-linear conjugate gradient initializer.
 *
 * This function initializes the nlcg_ws object. As long as the
 * argument `n` is less than the max memory size allocated
 * `g->max_size`, then it simply sets the parameters for the objective
 * function. It returns a success or failure code.
 */
int nlcg_set( objective_fn f, unsigned int n, void *p, nlcg_ws g){
  if( n > g->max_size ){
    return NLCG_MEM_ERROR;
  }
  g->lf.of.f = f;
  g->lf.of.n = n;
  g->lf.of.p = p;
  return NLCG_SUCCESS;
}

/*!
 * Non-linear conjugate gradient destructor.
 *
 * This function frees the memory allocated in the `nlcg_malloc`
 * function.
 */
void nlcg_free( nlcg_ws g){
  free(g->lf.dfdx);
  free(g->lf.s);
  free(g);
}
  

/*!
 * Non-linear conjugate gradient optimization.
 *
 * This function optimizes the objective function that is contained in
 * the workspace object `g`. The inputs are the vector (array) `x` and
 * the workspace `g`. The algorithm uses the vector `x` as a starting
 * point for the algorithm. 
 *
 * After optimization, the function returns the value of the objective
 * function at the minimum, and the vector `x` is at the position of
 * the minimum.
 */
double nlcg_optimize( double *x, nlcg_ws g){
  int i, j;
  double f, denom, num, beta, slope;
  int n = g->lf.of.n;
  opt_fn of = g->lf.of;
  lin_fn lf = g->lf;
  double *dfdx = lf.dfdx;
  double *s = lf.s;
  double *dfdx_old = g->dfdx_old;
  lf.x = x;

  /* Start off at the initial position */
  f = opt_fn_eval( x, dfdx, &of);
  denom = 0.;
  for( i=0; i<n; i++){
    denom += dfdx[i]*dfdx[i];
    /* The initial search direction is the negative gradient */
    s[i] = -dfdx[i];
  }
  for( j=0; j<NLCG_ITER_MAX; j++ ){
    /* Do a line search*/
    memcpy( dfdx_old, dfdx, n*sizeof(double));
    f = sw_line_search( f, &lf);
    /* Calculate the square magnitude of the gradient */
    printf( "value: %1.3e \n", f);
    slope = dfdx[0]*dfdx[0];
    num = dfdx[0]*(dfdx[0]-dfdx_old[0]);
    for( i=1; i<n; i++){
      slope += dfdx[i]*dfdx[i];
      num += dfdx[i]*(dfdx[i]-dfdx_old[i]);
    }
    /* Stop if falls below some tolerance */
    if( slope < DFDX_TOL*DFDX_TOL ){
      break;
    }
    /* Calculate the beta factor */
    beta = num / denom;
    beta = beta > 0 ? beta : 0;
    denom = slope;
    slope = 0.;
    /* Use beta factor to get new search direction */
    for( i=0; i<n; i++){
      s[i] = -dfdx[i] + beta*s[i];
      slope = s[i]*dfdx[i];
    }
    if( slope > 0. ){
      for( i=0; i<n; i++){
	s[i] = -dfdx[i];
      }
    }
    /* And repeat */
  }

  return f;
}

/*!
 * Line search function.
 *
 * This function performs the line search of the objective function
 * from a starting point `x0` in the search direction `s`. It searches
 * until the function is _sufficiently minimized_, where in this case
 * it means that it satisfies the strong Wolfe conditions. The line
 * search should need as inputs

 * - `f0` the value of the function at the starting point,
 * - `df0` the derivative in the search direction at the starting point,
 * - `phi(a)` the linear function related to the objective function.

 * The value at the starting point `f0` is given by the argument
 * `f`. The derivative `df0` can be calculated from the gradient
 * `dfdx` and the search direciton `s`, both of which are contained in
 * the linear function helper struct `lf` passed as the second
 * argument. In addition, `lf` contains all the information neccesary
 * for the linear function `phi(a)`, including the objective funciton
 * in the optimization function helper struct `of`, the starting point
 * `x`, and the search direction `s`.
 *
 * After _sufficiently minimizing_, the line search algorithm returns
 * the value of the function at the minimized point and
 * the vectors `x` and `dfdx` in the ls function contain the point at
 * which it is minimized, and the gradient at that minimum.
 */
static double sw_line_search( double f, lin_fn *lf){
  int i;
  double c1=C1, c2=C2;
  double f_0, df_0, a_m, f_m, df_m, a_p, f_p, df_p;
  int n = lf->of.n;
  double *s = lf->s;
  double *dfdx = lf->dfdx;
 
  /* set the initial function value */
  f_0 = f;
  /* set the initial value of a */
  lf->a_prev = 0.;
  /* cacluate the derivative and the magitude of the `s` vector */
  lf->smag = s[0]*s[0];
  df_0 = dfdx[0]*s[0];
  for( i=1; i<n; i++){
    lf->smag += s[i]*s[i];
    df_0 = dfdx[i]*s[i];
  }
  lf->smag = sqrt( lf->smag);
  df_0 /= lf->smag;

  /* set the smaller value `a_m` to zero */
  a_m = 0.;
  /* set the function and derivative at the smaller `a_m` */
  f_m = f_0;
  df_m = df_0;
  /* set the initial step length */
  a_p = ALPHA_0;
  for( i=0; i<LS_ITER_MAX; i++){
    /* evaluate the next value */
    f_p = lin_fn_eval( a_p, &df_p, lf);
    /* check if bracketed & `f_p` > `f_m` */
    if( f_p > f_0 + c1*a_p*df_0 || (f_p > f_m && i > 0)){
      return sw_bracket_search( f_0, df_0, a_m, f_m, df_m,
				a_p, f_p, df_p, lf);
    }
    /* check if satisfies the strong Wolfe condition */
    if( fabs(df_p) <= -c2*df_0 ){
      return f;
    }
    /* check if bracketed & `f_m` > `f_p` */
    if( df_p >= 0. ){
      return sw_bracket_search( f_0, df_0, a_p, f_p, df_p,
				a_m, f_m, df_m, lf);
    }
    /* move to larger range */
    a_m = a_p;
    f_m = f_p;
    df_m = df_p;
    a_p *= 2.;
  }
  /* if it always fails return the last value found */
  return f;
}

/*!
 * Bracketed line search function.
 *
 * This function is the second part of the line search function
 * `sw_line_search`. This function is called when it is known that
 * there exist a value \f$ \alpha_{\text{lo}} < \alpha <
 * \alpha_{\text{hi}} \f$ that statisfies the strong Wolfe
 * conditions. It takes as arguments

 * - `f_0` the value \f$ \phi(\alpha=0\) \f$,
 * - `df_0` the value \f$ \phi'(\alpha=0\) \f$,
 * - `a_lo` the bracketing value \f$ \alpha_{\text{lo}} \f$,
 * - `f_lo` the value \f$ \phi(\alpha_{\text{lo}}\) \f$,
 * - `df_lo` the value \f$ \phi'(\alpha_{\text{lo}}\) \f$,
 * - `a_hi` the bracketing value \f$ \alpha_{\text{hi}} \f$,
 * - `f_hi` the value \f$ \phi(\alpha_{\text{hi}}\) \f$,
 * - `df_hi` the value \f$ \phi'(\alpha_{\text{hi}}\) \f$,
 * - `lf` the linear function helper struct.

 * Identical to the `sw_line_search`, the bracket search algorithm
 * returns the value of the function at the minimized point, as well
 * as the minizing point and gradient in the `lf` struct.
 */
static double sw_bracket_search( double f_0, double df_0,
				 double a_lo, double f_lo, double df_lo,
				 double a_hi, double f_hi, double df_hi,
				 lin_fn *lf){
  int i;
  double c1=C1, c2=C2;
  double a_j, f_j, df_j;
  for( i=0; i<BRACKET_ITER_MAX; i++){
    /* interpolate between the two points */
    a_j = 0.5*(a_lo + a_hi);
    /* evaluate at the new point */
    f_j = lin_fn_eval( a_j, &df_j, lf);
    /* if appropriate, replace hi point*/
    if( f_j > f_0 + c1*a_j*df_0 || f_j > f_lo ){
      a_hi = a_j;
      f_hi = f_j;
      df_hi = df_j;
    }else{
      /* if satisfies strong Wolfe conditions, return */
      if( fabs(df_j) <= -c2*df_0){
	return f_j;
      }
      /* if appropriate make the lo point the new hi point */
      if( df_j*(a_hi-a_lo) >= 0. ){
	a_hi = a_lo;
	f_hi = f_lo;
	df_hi = df_lo;
      }
      /* replace lo point. */
      a_lo = a_j;
      f_lo = f_j;
      df_lo = df_j;
    }
  }
  /* if that doesn't finish after BRACKET_ITER_MAX turns, */
  /* return the point that evaluated smallest evaluated value */
  if( f_j < f_lo ){
    return f_j;
  }else{
    return lin_fn_eval( a_lo, &df_j, lf);
  }
}

/*!
 * Optimization function evaluation.
 *
 * This function evaluates the objective function at the vector `x`,
 * and return the function value, and the gradient `dfdx`. It uses the
 * optimization function helper struct `of` to hold all the extraneous
 * stuff.
 */
static double opt_fn_eval( const double *x, double *dfdx, const opt_fn *of){
  return of->f( of->n, x, dfdx, of->p);
}

/*!
 * Linear function evaluation.
 *
 * This function evaluates the function,

 * \f[ \phi(\alpha) = f({\bf x} + alpha {\bf s}), \]

 * It returns the value of of the function, as well as the derivative
 * with respect to \f$ \alpha \f$. Note: the function treats it as
 * though the search direction `s` is always a unit vector.
 */
double lin_fn_eval( double a, double *dfda, lin_fn *lf){
  int i;
  double f;
  int n = lf->of.n;
  double *x = lf->x;
  double *dfdx = lf->dfdx;
  double *s = lf->s;
  /* normalise the step length as if `s` were a unit vector */
  double step = (a-lf->a_prev)/lf->smag;
  /* move `x` by the step length */
  for( i=0; i<n; i++){
    x[i] = x[i] + step*s[i];
  }
  /* evaluate function */
  f = opt_fn_eval( x, dfdx, &(lf->of));
  /* calculate new dervative */
  *dfda = dfdx[0]*s[0];
  for( i=1; i<n; i++){
    *dfda += dfdx[i]*s[i];
  }
  *dfda /= lf->smag;
  /* set the previous `a` value */
  lf->a_prev = a;
  /* return function value */
  return f;
}
