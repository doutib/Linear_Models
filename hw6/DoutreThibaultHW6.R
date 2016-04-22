
## # Load data
data = read.table('mle.txt')
data = unlist(data, use.names = FALSE)
n = length(data)

## # Compute log likelihood
log_likelihood_aux = function(theta){
  n*log(theta) - 2*sum(log(theta+data))
}
log_likelihood = Vectorize(log_likelihood_aux)

# Plot
thetas = seq(from = 0.1, to = 100, by = 0.1)
lls = log_likelihood(thetas)
plot(thetas,lls,
     main = "Log likelihood of the data",
     ylab = "Log Likelihood",
     xlab = "Theta",
     type="l")

## # Change variable theta<-exp(phi)
log_likelihood_phi_aux = function(phi){
  n*phi - 2*sum(log(exp(phi)+data))
}
log_likelihood_phi = Vectorize(log_likelihood_phi_aux)

# Plot
phis = seq(from = -10, to = 10, by = 0.1)
lls = log_likelihood_phi(phis)
plot(phis,lls,
     main = "Log likelihood of the data",
     ylab = "Log Likelihood",
     xlab = "Phi",
     type="l")

## # Compute first and second derivatives
D_log_likelihood_phi_aux = function(phi){
  n - 2*sum(1/(exp(-phi)*data+1))
}
D_log_likelihood_phi = Vectorize(D_log_likelihood_phi_aux)

D2_log_likelihood_phi_aux = function(phi){
  - 2*exp(-phi)*sum(data/(exp(-phi)*data+1)^2)
}
D2_log_likelihood_phi = Vectorize(D2_log_likelihood_phi_aux)


## # Optimization without derivatives

# Function to optimize
f = function(phi){
  -log_likelihood_phi_aux(phi)
}

opt = optim(0,f,
            method="Brent",
            lower = 0,
            upper = 10)

points(opt$par,
       -opt$value,
       col = "red")

## # Optimization using derivatives

# # Compute first derivative
D_f = function(phi){
  -n + 2*sum(1/(exp(-phi)*data+1))
}
# # Compute second derivative
D2_f = function(phi){
  + 2*exp(-phi)*sum(data/(exp(-phi)*data+1)^2)
}

# # Newton raphson: one variable
newton_raphson = function(start,f,Df,D2f,eps){
  x_old=start-2*eps
  x_new=start
  while (abs(x_new-x_old)>eps){
    x_old = x_new
    x_new = x_old - Df(x_old)/D2f(x_old)
  }
  return(list(par = x_new, value = f(x_new)))
}

eps=10e-5
nr = newton_raphson(0,f,D_f,D2_f,eps)

points(nr$par,
       -nr$value,
       col = "green",cex=.8,pch=20)

theta_hat = exp(nr$par)
theta_hat

## # Standar error on theta

D2_log_likelihood_theta = function(theta){
  -n/theta^2+2*sum(1/(theta+data)^2)
}

# # Variance
var_hat = -1/D2_log_likelihood_theta(theta_hat)
# # Standard Error
se_hat = sqrt(var_hat)
se_hat

