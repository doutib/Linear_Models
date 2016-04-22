########### ### PART 1

## # Paramaters
n = 36
p = 3

## # Compute A = X'X
A = matrix(c(1,0.52,0.52,1),ncol=2)*n
beta_hat = solve(A,c(-0.26,-0.42)*n)
beta_hat

## # Compute adjusted variance
var_hat = n/(n-p)*(1-beta_hat[1]^2-beta_hat[2]^2-
                     2*beta_hat[1]*beta_hat[2]*0.52)
var_hat


## # Compute adjusted SE
# SE for both beta coefficients 
SE = sqrt(var_hat*solve(A)[1,1])
SE

## # Compute adjusted SE
V = var_hat*solve(A)
var_diff = V[1,1]+V[2,2]-2*V[1,2]
SE_diff = sqrt(var_diff)
diff = beta_hat[1]-beta_hat[2]
SE_diff


########### ### PART 2


## # Student test for beta 1.
SE1 = sqrt(var_hat*solve(A)[1,1])
t1 = beta_hat[1]/SE1
2*(1-pt(abs(t1),df=n-p)) #pvalue

## # Student test for beta 2.
SE2 = sqrt(var_hat*solve(A)[2,2])
t2 = beta_hat[2]/SE2
2*(1-pt(abs(t2),df=n-p)) #pvalue

## # Student test for the difference of beta 1 and 
t = (beta_hat[1]-beta_hat[2])/SE_diff
1-pt(abs(t),df=n-p) #pvalue



###### ## Divers (useless here)

## # Generate correlated data
corr_data = function(numobs, R, names=c()){
  # numobs: number of observations to generate
  # R: correlation matrix (square)
  U = t(chol(R))
  nvars = dim(U)[1]
  numobs = 100000
  random.normal = matrix(rnorm(nvars*numobs,0,1), nrow=nvars, ncol=numobs);
  X = U %*% random.normal
  newX = t(X)
  data = as.data.frame(newX)
  if (length(names)== nrow(R)){
    names(data) = names
  }
  data
}