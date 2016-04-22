
### ## ### LAB 11 ### ## ###

# Compute log likelihood
log_likelihood_aux = function(theta,X){
  n = lengthgth(X)
  n*log(theta) - 2*sum(log(theta+X))
}
log_likelihood = Vectorize(log_likelihood_aux)

# Change variable theta<-exp(phi)
log_likelihood_phi_aux = function(phi,X){
  n = length(X)
  n*phi - 2*sum(log(exp(phi)+X))
}

MLE = function (X){
  f = function(phi){-log_likelihood_phi_aux(phi,X)}
  opt = optim(0,f,
              method="Brent",
              lower = 0,
              upper = 10)
  theta_hat = exp(opt$par)
  theta_hat
}

## # 1-3. Generate uniform RV

set.seed(1)
theta_hat_values = c()
for (i in 1:1000){
  # Generate data
  U = runif(50)
  theta = 25
  X = theta*sapply(U,function(x) x/(1-x))
  # Compute MLE
  theta_hat = MLE(X)
  theta_hat_values = c(theta_hat_values,theta_hat)
}

## # 4. Plot histogram
hist(theta_hat_values,
     xlab = "theta_hat",
     main = "Realizations of MLE")
  
## # 5. Mean/SD
mu = mean(theta_hat_values)
mu
sigma = sd(theta_hat_values)
sigma
  
fisher_info = function(theta){
  1/(3*theta^2)
}

# Comparison
asympt_var = 1/sqrt(50*fisher_info(25))
sigma
asympt_var-sigma

# The asymptotic sd is equal to 1/sqrt(...)
# Therefore they shoud be equal for an infinite 
# number of simulations
# c.f. formula Example 4, Chapter 7

## # 6. 
# The asymptotic sd should be more accurate to estimate SE
# because it is a limit of the sd of a n-sized sample.

## # 7.
# The asymptotic sd doubles and the fisher info is divided
# by 4, see formula. The standard deviation of the observed info
# will also double and will still tend to the asymptotic sd
# as n grows to infinity.

### ## ### LAB 12 ### ## ###

data = read.table("pac01.dat")
head(data)
names(data)=c("AGE",
              "SEX",
              "RACE",
              "ETHNICITY",
              "MARITAL",
              "NUMKIDS",
              "FAMPERS",
              "EDLEVEL",
              "LABSTAT",
              "CLASSWORK",
              "FULLPART",
              "HOURS",
              "WHYNOTWORK",
              "INSCHOOL",
              "INDUSTRY",
              "OCCUPATION",
              "PINCOME",
              "INCFAM",
              "CITIZEN",
              "IMMIGYR",
              "HHSEQNUM",
              "FSEQNUM",
              "PERSCODE",
              "SPOUCODE",
              "FINALWGT",
              "MARCHWGT",
              "STATE")

head(data)
subset = data[,c("LABSTAT","AGE","SEX","RACE","EDLEVEL")]
head(subset)

# Split age
split_age1 = rep(0,length(data$AGE))
split_age2 = rep(0,length(data$AGE))
split_age3 = rep(0,length(data$AGE))

split_age1[data$AGE>=20 & data$AGE<=39] = 1#"20–39"
split_age2[data$AGE>=40 & data$AGE<=64] = 1#"40-64"
split_age3[data$AGE>=65] = 1#"65+"

# Split Race
split_race = rep(0,length(data$RACE))
split_race[data$RACE!=1] = 1#"white"

# Split Education level
split_edlevel1 = rep(0,length(data$EDLEVEL))
split_edlevel2 = rep(0,length(data$EDLEVEL))
split_edlevel1[data$EDLEVEL==39] = 1#"HS "
split_edlevel2[data$EDLEVEL>=40] = 1#"HS+"

# Split SEX
split_sex = rep(NA,length(data$SEX))
split_sex[data$SEX==1]=0#"M"
split_sex[data$SEX==2]=1#"F"

features = data.frame(LABSTAT = as.numeric(data$LABSTAT==1),
                      SEX = as.numeric(split_sex),
                      AGE1 = as.numeric(split_age1),
                      AGE2 = as.numeric(split_age2),
                      AGE3 = as.numeric(split_age3),
                      RACE = as.numeric(split_race),
                      EDLEVEL1 = as.numeric(split_edlevel1),
                      EDLEVEL2 = as.numeric(split_edlevel2))
head(features)
any(is.na(features))

## # 1. 
design = features[,-1]
design$Intercept = as.numeric(1)
head(design)
names_features = names(design)
design = as.matrix(design)


# Size of design matrix
dim(design)
Y=features[,1]

## # 2.
library(pracma)

Eta = function(x){
  sapply(x,function(x) 1/(1+exp(-x)))
}

# Negative log likelihood
loss = function(beta){
  x = design %*% beta
  -(sum(Y*log(Eta(x))+(1-Y)*log(1-Eta(x))))
}

beta0 = as.matrix(rep(0,dim(design)[2]))
loss(beta0)
mini = fminsearch(loss,beta0)
beta_hat = mini$xval
names(beta_hat) <- names_features
beta_hat

## Alternative
glm.fit <- glm(LABSTAT ~. , data = features, family = "binomial")
summary(glm.fit)

## # 3.
# Standard errors
summary(glm.fit)$coefficients[,2]

## # 4.
# When looking at the sign of the coefficients we can first
# say that employment is positively correlated with having
# been to high school or above, since the baseline is no
# high school and the coefficients for EDLEVEL1 and 
# EDLEVEL2 are positive. Similarly, we can also conclude 
# that beign either a woman or non-white has a bad impact 
# on employment.
# Moreover, it is important to notice that the p-values are 
# all very small, which shows the importance of all the
# features used to predict the outcome.


## # 5.
# First of all there is no way to correctly quantify the
# education, so we cannot use a real valued variable and
# perform regression. As a matter of fact we have categorical
# variables. We use dummy variables in order to avoid giving
# more weight to some variables: the way of assigning factors
# matters in the regression in the sense that the linear model
# gives more importance to categories with a bigger factor.


## # 6.
# The fact that most of women give birth may impact their
# employment since the employers know that they might not be
# able to work for a while.
# Moreover, the SEX variable might be correlated with the
# edication level for example. It is known that women have
# less access to education that men for some reasons and it
# might impact their employment.
# LABSTAT codes more than 4 are relevant because women are 
# more likely not to work than men do. Therefore they might
# not be looking for a job, which is generally not the case
# for men.
  
  