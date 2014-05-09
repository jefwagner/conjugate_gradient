Non-linear conjugate gradient optimization
==========================================

About a year and a half ago, I started a research project that needed
to optimize a function of up to 9000 variables. I could calculate the
gradient of this function, but not the Hessian (the matrix of second
derivatives). With so many variables, I decided to use a simple
optimization routine and settled on either a Non-linear conjugate
gradient
[https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method] or
a limited memory BFGS method
[https://http://en.wikipedia.org/wiki/Limited-memory_BFGS].

There were two C packages that I found that implemented these
routines, 

  1. The Gnu Scientific Library [http://www.gnu.org/software/gsl/], 
  2. nlopt from the Ab-initio group at MIT 
     [http://ab-initio.mit.edu/wiki/index.php/NLopt].  

However in testing both, of these packages would fail without any
indication as to why they failed. So I decided to write my own
optimization routine. Because it was simpler, I decided to implement
the non-linear conjugate gradient. Please find the details below.

Object oriented style 
--------------------- 

In order to optimize a function of `N` variables the non-linear
conjugate gradient routine needs memory to hold two vectors of lenth
`N` to hold the gradient and a search direction.

I implemented the memory managment by including the two pointers in an
opaque pointer as below

````c
struct cg_workspace_struct{
  double *dx /* gradient */
  double *s /* search direction */
  /* other stuff */
};
typedef struct cg_workspace_struct* cg_workspace;
````

and created two functions that allocate the memory and free the memory
in the opaque pointer.  

````c 
cg_workspace cg_malloc( unsigned int n);
void cg_free( cg_workspace g); 
````

The optimization function
-------------------------

I want the routine to be 