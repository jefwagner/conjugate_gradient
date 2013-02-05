#include <stdio.h>
#include <stdlib.h>
#include "ls.h"
#include "testfunc.h"

static int num_evals;

double myfunc( int n, const double *x, double *dx, void *p)
{
  int i;
  num_evals++;


  fprintf( stdout, "*** ");
  for( i=0; i<n; i++){
    fprintf( stdout, "%1.3g ", 1-x[i] );
  }
  fprintf( stdout, "\n");

  return( rosenbrock( n, x, dx, p));
}

int main( int argc, char **argv){
  int i, n;
  struct cg_workspace *g;
  double *x;

  if( argc != 2 ){
    fprintf( stderr, "Incorrect arguments");
    return 1;
  }
  n = atoi( argv[1]);

  x = (double *) malloc( n*sizeof(double) );

  num_evals = 0;
  g = allocate_cg_workspace( n);
  g->h = &myfunc;
  g->params = NULL;

  for( i=0; i<n; i++){
    x[i] = 0.;
  }

  fprintf( stdout, " deltax dhdx \n" );
  
  conjugate_gradient( g, x);
  
  fprintf( stdout, "Program finished after %d evaluations\n", num_evals);
  fprintf( stdout, "Found final point at:\n");
  for( i=0; i<n; i++){
    fprintf( stdout, "%1.3g ", x[i] );
  }
  fprintf( stdout, "\n");


  free_cg_workspace( g);
  free( x);

  return 0;
}
  
  
