/*!
 * Conjugate Gradient data structure
 * contains:
 * + common elements:
 *   - function pointer
 *   - size of function
 *   - functin parameter
 * + memory for arrays
 *   - derivative for current point
 *   - search direction
 *   - derivative for previous point
 *   - derivative for next point
 */
struct cg_workspace{
  double (*h)( int n, const double *x, double *dx, void *p);
  int n;
  void *params;

  /* variables used in the conjugate gradient method */
  double h0; /* value at current point */
  double *dx0; /* current gradient */
  double *s; /* search direction */

  /* variables used in the line search */
  double smag;
  double dh0;  /* derivative in search direction at current point */
  double alpham; /* step length at step i-1 */
  double hm; /* value at step i-1 */
  double dhm; /* derivative in search direction at step i-1 */
  double alphap; /* step length at step i+1 */
  double hp; /* value at step i */
  double *dxp; /* gradient at step i */
  double dhp; /* derivative in search direction at step i */
};

struct cg_workspace* allocate_cg_workspace( int n);
void free_cg_workspace( struct cg_workspace *g);
void conjugate_gradient( struct cg_workspace *g, double *x);

