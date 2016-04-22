library(pls)
rm(list = ls())
cat("\014")

# Data --------------------------------------------------------------------

makedata=function(p=20,wh=15,n=100){
  # This function generates p independent X variables and then uses wh of them
  # to generate Y based on coefficients determined by exp(exps) below.
  # Outside of the function, the data frame named data is further changed to
  # randomize the order of the variables so it's hard to tell for sure which
  # ones were used to generate Y.  You can always go back and examine switch
  # to figure out if your models include the right variables and which ones.
  X=matrix(rnorm(n*p),n,p)
  exps=seq(-1,-2.5,length=wh)
  beta=rep(0,p)
  beta[1:wh]=exp(exps)
  Y=.5+X%*%beta+rnorm(n)
  data = data.frame(Y,X)
  colnames(data)=c("Y",letters[1:20])
  switch=sample(20)+1
  data=data[,c(1,switch)]
  names(beta)=names(data)[switch]
  return(list(df=data,beta=beta))
}

set.seed(1)
data=makedata()
data$beta

# Replicate ---------------------------------------------------------------
cv_pcr = function(n){
  data = makedata(n=n)
  df = data$df
  beta0 = data$beta
  pcr.fit = pcr(Y~0+., ncomp=20, data=df, validation ="CV")
  cv_error = as.data.frame(RMSEP(pcr.fit)$val)[1,]
  names(cv_error) = 0:20
  # beta pcr
  beta0 = data$beta
  beta_pcr = rev(sort(pcr.fit$coefficients[,,20]))
  return(list(CV=cv_error, beta0=beta0, beta=beta_pcr))
}


library(parallel)

pcr_rep = function(B,n=100){
  ncores = detectCores()-1
  print(paste('Starting ', ncores, ' cores...'))
  cl = makeCluster(ncores)
  clusterEvalQ(cl,library(pls))
  clusterExport(cl,list("cv_pcr","makedata"))
  R = parSapply(cl, 1:B, function(i,...) 
  { cv_pcr(n=n) } )
  CV = matrix(unlist(R["CV",]), nrow = B, byrow = T)
  beta0 = matrix(unlist(R["beta0",]), nrow = B, byrow = T)[1,]
  beta = matrix(unlist(R["beta",]), nrow = B, byrow = T)
  stopCluster(cl)
  print('Done.')
  return(list(CV=CV,beta=beta,beta0=beta0))
}

B = 100 # TODO: Change to 1000
Rep_pcr100 = pcr_rep(B,n=100)
Rep_pcr1000 = pcr_rep(B,n=1000)


# Plot --------------------------------------------------------------------

# n=100
boxplot(Rep_pcr100$CV, main=paste("CV error for PCR with n = ",100),
        xlab = "Number of components", ylab="test error",names=0:20)

#n=1000
boxplot(Rep_pcr1000$CV, main=paste("CV error for PCR with n = ",1000),
        xlab = "Number of components", ylab="test error",names=0:20)

#n=100
boxplot(Rep_pcr100$beta, main=paste("beta_hat for PCR with n = ",100),
        xlab = "Number of components", ylab="beta")
points(1:20,Rep_pcr1000$beta0,col="red")

#n=1000
boxplot(Rep_pcr1000$beta, main=paste("beta_hat for PCR with n = ",1000),
        xlab = "Number of components", ylab="beta")
points(1:20,Rep_pcr1000$beta0,col="red")


