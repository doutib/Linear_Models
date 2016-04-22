###  HW 3   Doutre Thibault
rm(list = ls())

# Question 1 --------------------------------------------------------------

generate_data = function(n = 100){
  # Define parameters
  rho = 0.7
  mu1=180; s1=40; mu2=66; s2=3
  
  # Define X, Y and Z with the bivariate normal relationship
  X = rnorm(n)
  Z = rnorm(n)
  eps = sqrt(1-rho^2) * Z
  Y = rho * X + eps
  
  # Adjust means and variances
  Y = (Y-mean(Y))/sd(Y)*s2+mu2
  X = (X-mean(X))/sd(X)*s1+mu1
  
  # Adjust rho by transforming Y
  rho_hat = cor(X,Y)
  a = s1^4*(rho^2-1)
  b = 2*rho_hat*s1^3*s2*(rho^2-1)
  c = (rho^2-rho_hat^2)*s2^2*s1^2
  delta = b^2-4*a*c
  correction = (-b-sqrt(delta))/(2*a)
  Y=Y+correction*X
  
  # Adjust means and variances
  Y = (Y-mean(Y))/sd(Y)*s2+mu2
  X = (X-mean(X))/sd(X)*s1+mu1
  
  # Put into data frame
  df = data.frame(WT = Y,
                  HT = X,
                  BMI = 703 * Y / X^2)
  
  # Output
  return(list(df=df,rho=rho,eps=eps))
}
data = generate_data()

df = data$df
M=df[,1:2]
rho = data$rho
eps = data$eps

# Correlation and covariance matrices
cor(M) 
cov(M)

# Mean of variables
mean(df$WT)
mean(df$HT)


# Question 2 --------------------------------------------------------------

lm.fit = lm(WT ~ ., data = df)
beta = lm.fit$coefficients

par(mfrow=c(2,2))
plot(lm.fit)
par(mfrow=c(1,1))

# In particular, the sd of the residuals is not equal to 1:
sd(lm.fit$residuals)


# Question 3 --------------------------------------------------------------


# True value of beta
beta_true = rho*sd(df$WT)/sd(df$HT)
beta_true

# The simulated value of $\beta_1$ is
beta1 = beta["HT"]
beta1


# Question 4 --------------------------------------------------------------


# Question 5 --------------------------------------------------------------

WT_hat = predict(lm.fit)
e = df$WT-WT_hat

# To see if the variables are correlated, wee look for some pattern in the plot of one against the other.
plot(eps,df$HT) # not correlated, seem to be independent
plot(eps,df$WT) # positively correlated -> dependent
plot(e,df$WT) # positively correlated -> dependent
plot(e,df$HT) # not correlated, seem to be independent

# To see if two vectors are orthogonal, I compute their scalar product.
sum(e*df$WT) # not orthogonal
sum(e*df$HT) # orthogonal
sum(eps*df$WT) # not orthogonal
sum(eps*df$HT) # not orthogonal

# Question 6 --------------------------------------------------------------

# Generate beta
beta = function(n=100){
  data0 = generate_data(n)
  df0 = data0$df
  lm.fit0 = lm(WT ~ HT+BMI, data = df0)
  lm.fit0$coefficients
}

# Replication
betas = replicate(1000,beta()["HT"])

#  Question 7 -------------------------------------------------------------

hist(betas)
abline(v = beta_true) 
# bias
mean(betas) - beta_true
