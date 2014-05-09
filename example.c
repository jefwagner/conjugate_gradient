/*!
 * Example file on how to use the Conjugate gradient function
 * contained in the file cg.c and cg.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include "cg.h" /*! include the header file */


/*!
 * This is a standard test function called the rosenbrock function. It
 * has a narrow curved vally in n dimensional space. The minimum is
 * located in that vally at \f[x_i=1 \forall i\f]. This function is a
 * fairly standard function for testing optimization routines. The
 * parameter that is passed simply holds an integer that is increased
 * by 1 every time the function is evaluated.
 */
double rosenbrock( int n, const double *x, double *dx, void *p)
{
  int i;
  double sum=0.;

  /*! increment the integer held in p by 1 */
  int *count = (int *) p;
  (*count)++;

/*   for( i=0; i<n; i++){ */
/*     fprintf( stdout, "%1.3g ", x[i] ); */
/*   } */
/*   fprintf( stdout, "\n"); */

  sum += (1.-x[0])*(1.-x[0]) +
    100.*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
  dx[0] = -2.*(1-x[0]) -
    400*x[0]*(x[1]-x[0]*x[0]);
  for( i=1; i<(n-1); i++){
    sum += (1.-x[i])*(1-x[i]) + 
      100.*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]);
    dx[i] = -2.*(1.-x[i]) -
      400.*x[i]*(x[i+1]-x[i]*x[i]) +
      200.*(x[i]-x[i-1]*x[i-1]);
  }
  dx[n-1] = 200.*(x[n-1]-x[n-2]*x[n-2]);

  return( sum);
}

/*!
 * The example program shows how to use the conjugate gradient
 * algroithm. It takes as argument a single integer that determines
 * the dimension of the system to be minimized is.
 *
 * The program has the following steps
 * 1. Reads the argument from the command line
 * 2. Allocates the memory needed
 * 3. Assigns the function and parameter to the minization object
 * 4. Initialized the first guess
 * 5. Finds the minimum using the Nonlinear CG algroithm
 * 6. Output the results
 * 7. Frees up the memeory allocated
 */
int main( int argc, char **argv){
  int i, n, count;
  struct cg_workspace *g;
  double *x;

  if( argc != 2 ){
    fprintf( stderr, "Incorrect arguments");
    return 1;
  }
  n = atoi( argv[1]);

  x = (double *) malloc( n*sizeof(double) );
  g = allocate_cg_workspace( n);

  g->h = &rosenbrock;
  count = 0;
  g->params = (void *) &count;

  for( i=0; i<n; i++){
    x[i] = 0.;
  }

  conjugate_gradient( g, x);
  
  fprintf( stdout, "Program finished after %d evaluations\n", count);
  fprintf( stdout, "cg counted  %d evaluations\n", g->eval);
  fprintf( stdout, "Found final point at:\n");
  for( i=0; i<n; i++){
    fprintf( stdout, "%1.3g ", x[i] );
  }
  fprintf( stdout, "\n");

  free_cg_workspace( g);
  free( x);

  return 0;
}
