Non-linear conjugate gradient optimization
==========================================

An early problem in calculas is to find a minimum or maximum value of
a (often 1-dimentional) function.

We are attempting to find a minimum or maximum value of and
n-dimentional function.

Need a function that returns value f and gradient df at a point x.
(f, df) = f( x)

Non-linear conjugate gradient 
Find a point x_min where f(x) is near a minimum
Start at a point x
(f, df) = f( x)
s = -df
beta = df.df
(x, f, df) = min( f() , (x, s))
while x is not near min
    beta = min( 0, df.df/beta);
    s = df + beta*s
    (x, f, df) = min( f() , (x, s))
return x_min = x

Need a function that 

Line search
Find a value alpha such that f(x+alpha*s) is sufficiently minimized
(f_0, df_0) = f(x)
alpha_left = 0;
Choose alpha_right;
(f, df) = f(x+alpha_right*s)
loop
  if f >= f_0+c_1*alpha_right*ds
    return bracket_search( f(), (x, s), alpha_left, alpha_right);
  else if |df| <= -c_2*df0
    return (x+alpha*s, f, df)
  else if f < f_0
    return bracket_search( f(), (x, s), alpha_right, alpha_left);
  alpha_m = alpha_p
  alpha_p = 2 alpha_p

Bracket Search
alpha_j = 1/2( alpha_p+alpha_m)
(f, df) = f(x+alpha_j)

