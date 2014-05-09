/*!
 * Nonlinear Conjugate Gradient File
 */

struct nlcg_workspace_struct{
  double *dx;
  double *s;
};
typedef struct nlcg_workspace_struct* nlcg_workspace;

typedef double (*opt_func)( int n, const double *x, 
			    double *dfdx, void *params);

double nlcg_optimize( nlcg_workspace g, opt_func func,
		      unsigned int n, double *x, void *params){
  unsigned int i;
  double f, denom, beta, slope;
  double alpha = ALPHA_0;
  double *dx = g->dx;
  double *s = g->s;

  f = func( n, x, dx, params);
  denom = 0.;
  for( i=0; i<n; i++){
    denom += dx[i]*dx[i];
    s[i] = -dx[i];
  }

  f = sw_line_search( g, f, &alpha, func, n, x, params);
  slope = 0.;
  for( i=0; i<n; i++){
    slope += dx[i]*dx[i];
  }
  while( slope < DX_TOL*DX_TOL ){
    beta = slope / denom;
    beta = beta > 0 ? beta : 0;
    for( i=0; i<n; i++){
      s[i] = -dx[i] + beta*s[i];
    }
    f = sw_line_search( g, f, &alpha, func, n, x, params);
  }

  return f;
}

double sw_line_search( nlcg_workspace g, double f0,
		       double *alpha, opt_func func,
		       unsigned int n, double *x, void *params){
  double smag, dx0, alpham, fm, dxm, fp, dxp;
  double *dx = g->dx;
  double *s = g->s;

  smag = 0;
  for( i=0; i<n; i++){
    smag += s[i]*s[i];
  }
  smag = sqrt( smag);

  dx0 = 0.;
  for( i=0; i<n; i++){
    dx0 = dx[i]*s[i];
  }
  dx0 /= smag;
  alpham = 0.;
  fm = f0;
  dfm = dx0;

  while( *alpha < ALPHA_MAX ){
    for( i=0; i<n; i++){
      x[i] += *alpha*s[i];
    }
    fp = func( n, x, dx, params);
    dxp = 0.;
    for( i=0; i<n; i++){
      dxp = dx[i]*s[i];
    }
    dxp /= smag;
    
    if( fp > f0+c1*(*alpha)*dx0 || fp >= fm ){
      return ws_bracket_search( g, &alpha, 
				alpham, fm, dfm,
				*alpha, fp, dfp,
				func, n, x, params);
    }
    if( fabs( dxp) <= - c2*dx0 ){
      return fp;
    }
    if( dxp >= 0. ){
      return ws_bracket_search( g, &alpha,
				*alpha, fp, dfp,
				alpham, fm, dfm,
				func, n, x, params);
    }
    alpham = *alpha;
    fm = fp;
    dfm = dfp;
    *alpha *= 2.;
  }
  return fp;
}

double ws_bracket_search( nlcg_workspace g, double *alpha,
			  double a_lo, double f_lo, double df_lo,
			  double a_hi, double f_hi, doubel df_hi,
			  opt_func func, unsigned int n, 
			  double *x, void *params){
  *alpha = interpolate( a_lo, f_lo, df_lo, a_hi, f_hi, df_hi);
  for( i=0; i<n; i++){
    x[i] += *alpha*s[i];
  }
  f = func( n, x, dx, params);
  dx_mid = 0.;
  for( i=0; i<n; i++){
    dx_mid += dx[i]*s[i];
  }
  dx_mid /= smag;
  
  if( f > 
