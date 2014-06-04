/*!
 * Example file on how to use the Conjugate gradient function
 * contained in the file cg.c and cg.h.
 */
#include <stdio.h>
#include <stdlib.h>
/* include the header file */
#include "cg.h" 


/*!
 * Rosenbrock function. 
 *
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
 * Example nonlinear conjugate gradient optimization.
 *
 * The example program shows how to use the conjugate gradient
 * routine. The program has the following steps
 
 * 1. Reads the dimension from the command line
 * 2. Allocates the memory needed
 * 3. Assigns the function and parameter to the minization object
 * 4. Initialized the first guess
 * 5. Finds the minimum using the Nonlinear CG algroithm
 * 6. Output the results
 * 7. Frees up the memeory allocated
 */
int main( int argc, char **argv){
  int i, n, count;
  nlcg_ws g;
  double f, *x;

  /* completely useless error message*/
  if( argc != 2 ){
    fprintf( stderr, "Incorrect arguments!\n");
    return 1;
  }
  n = atoi( argv[1]);
  if( n < 2 ){
    fprintf( stderr, "Dimension must be 2 or larger.\n");
    return 1;
  }

  /* Allocate the memory */
  x = (double *) malloc( n*sizeof(double));
  g = nlcg_malloc( n);

  /* Populate the nlcg object */
  nlcg_set_sys( (objective_fn) &rosenbrock, n, (void *) &count, g);
  nlcg_set_tol( 0., 0., 1.e-3, 0, g);

  /* Randomize the initial guess */
  for( i=0; i<n; i++){
    x[i] = 4.*rand()/RAND_MAX -2.;;
  }

  /* Optimize the function */
  f = nlcg_optimize( x, g);
  
  /* Print the results */
  fprintf( stdout, "Optimization finished after %d evaluations\n", count);
  fprintf( stdout, "  Value at minimum : %1.6g\n", f);
  fprintf( stdout, "  Optimum point :\n   ");
  for( i=0; i<n; i++){
    fprintf( stdout, "%1.3g ", x[i] );
    if( (i+1)%6 == 0 ){
      fprintf( stdout, "\n   " );
    }
  }
  fprintf( stdout, "\n");

  /* Free up memory */
  nlcg_free( g);
  free( x);

  return 0;
}
