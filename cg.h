/* irene loves mommy  */
/* I love mommy and daddy! */

#ifndef JW_NLCG
#define JW_NLCG

#define NLCG_SUCCESS 1
#define NLCG_MEM_ERROR -1


/*!
 * typedef for an objective functions.
 *
 * This gives the prototype for the multidimensional objective function.
 * The input parameters:
 *
 * 1.) `n` : The dimension of the function,
 * 2.) `x` : An `n` dimension evaluation point,
 * 3.) `p` : A general pointer to point to other parameters.
 *
 * The output are:
 *
 * 1.) The return value : the value of the function at point `x`,
 * 2.) `dfdx` : An `n`, dimensional derivative 
 *
 * \f[ dfdx[i] = \frac{\partial f}{\partial x[i]}. \f]
 */
typedef double (*objective_fn)( int n, const double *x, 
				double *dfdx, void *p);

/*! nlcg object */
typedef struct nlcg_ws_struct* nlcg_ws;
/*! nlcg object constructor */
nlcg_ws nlcg_malloc( unsigned int max_size);
/*! destructor */
void nlcg_free( nlcg_ws g);
/*! tolerance setter */
void nlcg_set_tol( double df_tol, double dx_tol, double dfdx_tol, 
		   int max_eval, nlcg_ws g);
/*! system setter */
int nlcg_set_sys( objective_fn f, unsigned int n, void *p, nlcg_ws g);
/*! optimization routine */
double nlcg_optimize( double *x, nlcg_ws g);

#endif
