mle.dat = as.numeric(scan("mle.dat", what = "numeric"))
#### functions ####

# compute (log) likelihood
dty <- function( theta , x , log = FALSE ) {
  if( log ) {
    return( log(theta) - 2*log(theta+x) )
  } else {
    return( exp(log(theta) - 2*log(theta+x)) )
  }
}

# computes negative log likelihood
nll <- function( theta , x ) {
  return(-sum( dty( theta , x , log = TRUE ) ))  
}

# computes first order derivative of negative log likelihood
d1.nll <- function( theta , x ) {
  return(-(length(x)/theta - 2 * sum( 1/(theta+x) )))	
}

# computes second order derivative of negative log likelihood
d2.nll <- function( theta , x ) {
  return(-(-length(x)/theta^2 + 2 * sum( 1/(theta+x)^2 )))
}

# computes likelihood in term of phi
dty.phi <- function( phi , x , log = FALSE ) {
  return(dty( exp(phi) , x , log = log ))
}                                   

# computes negative log likelihood in term of phi
nll.phi <- function( phi , x ) {
  return(-sum( dty.phi( phi , x , log = TRUE ) ))
} 

# computes first order derivative of negative log likelihood in term of phi
d1.nll.phi <- function( phi , x ) {
  (-(length(x) - 2 * exp(phi) * sum( 1 / (exp(phi) + x) ) ))
}

# computes second order derivative of negative log likelihood in term of phi
d2.nll.phi <- function( phi , x ) {
  (-(-2 * exp(phi) * sum( x / (exp(phi) + x)^2 ) ))
}





#### optimization ####
optim( par = mean(mle.dat) , fn = nll.phi , 
       x = mle.dat , method = "BFGS" , 
       hessian = TRUE , control = list(trace=1) )

optim( par = log(mean(mle.dat)) , fn = nll.phi , 
	x = mle.dat , method = "BFGS" , 
	hessian = TRUE , control = list(trace=1) )

optim( par = log(mean(mle.dat)) , fn = nll.phi , 
	x = mle.dat , method = "L-BFGS-B" , 
	hessian = TRUE , control = list(trace=1) )

# all give the same result


# optimizing wrt theta
optim( par = mean(mle.dat), fn = nll, 
       x = mle.dat , method = "L-BFGS-B" , 
       hessian = TRUE , control = list(trace=1) )


# 