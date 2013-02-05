#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>

#define SQRT_2PI 2.50662827463

struct testfunc_params{
  double *sigma;
  double *x_avg;
};

double test_func( int n, const double *x, double *dx, void *params)
{
  struct testfunc_params * p = (struct testfunc_params *) params;
  double *sigma = p->sigma;
  double *x_avg = p->x_avg;
  int i;
  double sum_exp=0., prod_sig=1.;
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


  return( - f);
}

double gsl_f( const gsl_vector *x, void * params)
{
  double *dx = (double *) malloc( (x->size)*sizeof(double) );
  return test_func( x->size, x->data, dx, params);
  free (dx);
}

void gsl_df( const gsl_vector *x, void * params, gsl_vector *dx)
{
  test_func( x->size, x->data, dx->data, params);
}

void gsl_fdf( const gsl_vector *x, void *params,
	      double *f, gsl_vector *dx)
{
  *f = test_func( x->size, x->data, dx->data, params);
}

int main()
{
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  double x_avg[2] = {1.1534, -2.3441};
  double sigma[2] = {10, 15.25};

  struct testfunc_params p;
  p.x_avg = x_avg;
  p.sigma = sigma;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = 2;
  my_func.f = gsl_f;
  my_func.df = gsl_df;
  my_func.fdf = gsl_fdf;
  my_func.params = &p;

  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 0.0);
  gsl_vector_set (x, 1, 0.0);

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, 2);

  gsl_multimin_fdfminimizer_set ( s, &my_func, x, 0.01, 1.e-4);

  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status){ 
      fprintf( stderr, "error: %s\n", gsl_strerror (status) );
      break;
    }
    
    status = gsl_multimin_test_gradient (s->gradient, 1.e-10);

    if (status == GSL_SUCCESS){
      fprintf( stdout, "   ");
    }

    fprintf( stdout, "%1.3e %1.3e\n",
	     gsl_vector_get( s->x, 0),
	     gsl_vector_get( s->x, 1));
  }while( status == GSL_CONTINUE && iter < 100);

  gsl_multimin_fdfminimizer_free( s);
  gsl_vector_free( x);

  return 0;
}
