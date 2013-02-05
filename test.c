#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQRT_2PI 2.50662827463

struct cg_workspace{
  double (*h)( int n, const double *x, double *dx, void *p);
  int n;
  void *params;

  double val;
  double alpha;

  double *dx;
  double *dx_old;
  double *s;
  double *x_alpha;
  double *dx_alpha;

};

int conjugate_gradient( struct cg_workspace *g, double *x);
struct cg_workspace* 
allocate_cg_workspace( double (*h)(int n, const double *x, 
				   double *dx, void *p),
		       int n, void *params);
void free_cg_workspace( struct cg_workspace *g);

double test_func( int n, const double *x, double *dx, void *params);

struct testfunc_params{
  double *sigma;
  double *x_avg;
};

int main()
{
  struct testfunc_params p;
  double sigma[2], x_avg[2], x[2];
  int n=2;
  struct cg_workspace *cg;

/*   srand(0); */

/*   for( i=0; i<n; ++i){ */
/*     sigma[i] = 1.+10*( ((double) rand())/ ((double) RAND_MAX)); */
/*     fprintf( stdout, "%1.3e ", sigma[i]); */
/*     x_avg[i] = -10. + 20.*( ((double) rand())/ ((double) RAND_MAX)); */
/*     fprintf( stdout, "%1.3e\n", x_avg[i]); */
/*     x[i] = 0.; */
/*   } */

  sigma[0] = 10.; sigma[1] = 15.25;
  x_avg[0] = 1.1534; x_avg[1] = -2.3441;
  x[0] = 0.; x[1] = 0.;

  p.sigma = sigma;
  p.x_avg = x_avg;

  cg = allocate_cg_workspace( &test_func, n, (void *) &p);
  conjugate_gradient( cg, x);
  free_cg_workspace( cg);

  return( 0);
}

double test_func( int n, const double *x, double *dx, void *params)
{
  struct testfunc_params * p = (struct testfunc_params *) params;
  double *sigma = p->sigma;
  double *x_avg = p->x_avg;
  int i;
  double sum_exp=0, prod_sig=1;
  double f;

  for( i=0; i<n; ++i){
    prod_sig *= SQRT_2PI*sigma[i];
    sum_exp += (x[i]-x_avg[i])*(x[i]-x_avg[i])/(sigma[i]*sigma[i]);
  }
  prod_sig = 1./prod_sig;


  f = - 100.*prod_sig*exp( -0.5*sum_exp);
  for( i=0; i<n; ++i){
    dx[i] = -(x[i]-x_avg[i])/(sigma[i]*sigma[i])*f;
  }

  for( i=0; i<n; ++i){
    fprintf( stdout, "%1.3e ", x[i]);
  }
  for( i=0; i<n; ++i){
    fprintf( stdout, "%1.3e ", dx[i]);
  }
  fprintf( stdout, "\n");


  return( f);
}

