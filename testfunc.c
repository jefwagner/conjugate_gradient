#include <stdio.h>

double rosenbrock( int n, const double *x, double *dx, void *p)
{
  int i;
  double sum=0.;
  
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

double sphere( int n, const double *x, double *dx, void *p )
{
  int i;
  double sum = 0.;
  
  for( i=0; i<n; i++){
    sum += x[i]*x[i];
    dx[i] += 2.*x[i];
  }

  return sum;
}

double beale( int n, const double *x, double *dx, void *p )
{
  double c1, c2, c3;
  
  if( n != 2 ){
    fprintf( stderr, "Only defined for n=2");
  }
  
  c1 = 1.5-x[0]+x[0]*x[1];
  c2 = 2.25-x[0]+x[0]*x[1]*x[1];
  c3 = 2.625-x[0]+x[0]*x[1]*x[1]*x[1];

  return( c1*c1 + c2*c2 + c3*c3 );
}
