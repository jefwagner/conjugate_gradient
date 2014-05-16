/* irene loves mommy  */
/* I love mommy and daddy! */

#ifndef JW_NLCG
#define JW_NLCG

#define NLCG_SUCCESS 1
#define NLCG_MEM_ERROR -1

typedef double (*objective_fn)( int n, const double *x, 
				double *dfdx, void *p);
typedef struct nlcg_ws_struct* nlcg_ws;

nlcg_ws nlcg_malloc( unsigned int max_size);
int nlcg_set( objective_fn f, unsigned int n, void *p, nlcg_ws g);
double nlcg_optimize( double *x, nlcg_ws g);
void nlcg_free( nlcg_ws g);

#endif
