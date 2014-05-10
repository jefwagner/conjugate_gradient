/*!
 * Nonlinear Conjugate Gradient File
 */

typedef double (*opt_func)( int n, const double *x, 
			    double *dfdx, void *params);

struct nlcg_workspace_struct{
  unsigned int n;
  opt_func func;
  void *params;

  unsigned int array_size;
  double *dx; /*!< */
  double *s; /*!< */
};
typedef struct nlcg_workspace_struct* nlcg_workspace;

typedef struct{
  double smag;
  double alpha_m;
  double f_m;
  double df_m;
  double alpha_p;
  double f_p;
  double df_p;
} ls_workspace;

nlcg_workspace nlcg_malloc( unsigned int n){
  nlcg_workspace g = (nlcg_workspace) 
    malloc( sizeof(struct nlcg_workspace));
  g->array_size = n;
  g->dx = (double*) malloc( n*sizeof(double));
  g->s = (double*) malloc( n*sizeof(double));
}

int nlcg_set( nlcg_workspace g, unsigned int n, 
	      opt_func func, const void *params){
  if( n > g->array_size ){
    return NLCG_ARRAY_ERROR;
  }
  g->n = n;
  g->func = func;
  g->params = params;
  return NLCG_SUCCESS;
}

void nlcg_free( nlcg_workspace g){
  free( g->dx);
  free( g->s);
  free( g);
}

void ls_set( ls_workspace *ls,
	     double alpha_m, double f_m, double df_m,
	     double alpha_p, double f_p, double df_p){
  ls->alpha_m = alpha_m;
  ls->f_m = f_m;
  ls->df_m = df_m;
  ls->alpha_p = alpha_p;
  ls->f_p = f_p;
  ls->df_p = df_p;
}

double nlcg_optimize( nlcg_workspace g, double *x){
  unsigned int i;
  double f, denom, beta, slope;
  double *dx = g->dx;
  double *s = g->s;
  ls_workspace ls;

  f = func( n, x, dx, params);
  denom = 0.;
  for( i=0; i<n; i++){
    denom += dx[i]*dx[i];
    s[i] = -dx[i];
  }

  ls_set( &ls, 0., f, -sqrt( denom), 0., 0., 0.); 
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
