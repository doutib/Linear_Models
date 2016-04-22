### HW 4.  GLS and correlated errors.
### Part 2.  (Part 1 is the problems from the book.)
### 1) Generate correlated errors.

# First you write the code for making the matrix G.
# Use multiple lines if that seems easier.
# Then get the cholesky decomposition using the function chol().

G=
R=

# The next line makes everyone start at the same seed to get the same 
# random numbers.  For fun you may want to comment the line out to see 
# what happens for different simulations.

set.seed(12345)

# The following line will generate uncorrelated random errors.

z = matrix(rnorm(100000),100,1000)

# Using a loop may not be the fastest way to do this (or maybe it is), but 
# I think it is the easiest to understand and my goal is to make this 
# understandable.  It runs plenty fast in any case...

eps=matrix(0,100,1000)
for(i in 1:1000){
  eps[,i] = z[,i] %*% R
}

# Finish what you need to do in 1).

### 2. OLS

# Generate X

X = matrix(rnorm(100000),100,1000)

# That gets us past all the random number generation you'll need, I just
# wanted to make sure everyone would do this part in the same order
# so that everyone would get the same results.

# Now construct Y according to the model Y = Xb + eps
# where b=1 and Y will be a 100x1000 matrix, one column corresponding to
# each replication.  Note that when the HW refers to beta I'm talking about
# the coefficient b here.  Don't include an intercept term when you generate
# the data but do include one when you fit the model.

# Now you're on your own for the rest.  Don't use any fancy functions 
# (unless you want to check your work), I'd like you to do all this based 
# on matrix algebra for now.
