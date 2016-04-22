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
  return(data)
}


set.seed(1)
data=makedata()
lmout=lm(Y~.,data=data)
summary(lmout)
coef(summary(lmout))[-1,]



# Criteria ----------------------------------------------------------------

## # Cross validated MSE

MSE_cv = function(data, nfold=10, formula="Y~."){
  Y_cv = c()
  nrows = nrow(data)
  MSE_test = c()
  for (i in seq(0, nrows-nfold, by=nrows/nfold)){
    test = (i+1):(nfold+i)
    train = -test
    lm.fit = lm(formula, data=data[train,])
    Ytest = predict(lm.fit, data[test,])
    MSE_test = c(MSE_test, mean((data$Y[test]-Ytest)^2))
  }
  return(mean(MSE_test))
}

MSE_cv(data)


# Backward selection ------------------------------------------------------

backward_lm = function(data,nfold=10){
  # Initialize with OLS
  formula = "Y~."
  lm.fit = lm(formula, data=data)
  # Initialize outputs
  MSE_test = c()
  AIC = c()
  BIC = c()
  Cp = c()
  next_to_remove = ""
  variables = c()
  n = nrow(data)
  while(length(names(lm.fit$model))>1){
    MSE = MSE_cv(data, nfold=nfold, formula)
    MSE_test = c(MSE,MSE_test)
    AIC = c(AIC(lm.fit),AIC)
    BIC = c(BIC(lm.fit),BIC)
    d = length(names(lm.fit$coefficients))-1
    Cp = c(MSE + 2*d*mean((lm.fit$residuals)^2)/n,Cp)
    t_values = coef(summary(lm.fit))[,"t value"]
    # Variable with smallest t-value
    next_to_remove = names(which.min(t_values[-1]))
    # Store removed variables in the order
    variables = c(next_to_remove,variables)
    # Update formula
    formula = paste(formula,"-",next_to_remove,sep="")
    # Update model using new formula
    lm.fit = update(lm.fit, formula)
  }
  # Intercept only
  variables = variables
  MSE = MSE_cv(data, nfold, formula)
  MSE_test = c(MSE,MSE_test)
  AIC = c(AIC(lm.fit),AIC)
  BIC = c(BIC(lm.fit),BIC)
  Cp = c(MSE,Cp)
  return(list(variables = variables, MSE_test = MSE_test, 
              AIC=AIC, BIC=BIC, Cp=Cp))
}

backward = backward_lm(data)
backward

# Generate data and performs backward selection
backward_data = function(n=100,nfold=10){
  data = makedata(n=n)
  back=backward_lm(data,nfold=nfold)
  back$names=names(data)[-1]
  return(back)
}


## # REPLICATE
library(parallel)
metric_values = function(B,n=100){
  ncores = detectCores()-1
  print(paste('Starting ', ncores, ' cores...'))
  cl = makeCluster(ncores)
  clusterExport(cl,list("makedata", "backward_lm", "MSE_cv","backward_data"))
  R = parSapply(cl, 1:B, function(i,...) 
    { backward_data(n=n) } )
  MSE_values = matrix(unlist(R["MSE_test",]), nrow = B, byrow = T)
  AIC = matrix(unlist(R["AIC",]), nrow = B, byrow = T)
  BIC = matrix(unlist(R["BIC",]), nrow = B, byrow = T)
  Cp_values  = matrix(unlist(R["Cp",]), nrow = B, byrow = T)
  variables = matrix(unlist(R["variables",]), nrow = B, byrow = T)
  names = matrix(unlist(R["names",]), nrow = B, byrow = T)
  stopCluster(cl)
  print('Done.')
  return(list(MSE=MSE_values, AIC=AIC, BIC=BIC, Cp=Cp_values,
              variables=variables, names=names))
}

B = 100 # TODO: Change to 1000
Rep_backward100 = metric_values(B,n=100)
Rep_backward1000 = metric_values(B,n=1000)


# Plots -------------------------------------------------------------------

plot_err = function(metrics,n){
  m = sapply(metrics, colMeans)
  scale_m = apply(m, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
  plot(scale_m[,1],type="b",xlab="p",ylab="err",
       main=paste("Adjusted errors for n =",n), col="red")
  lines(type="b",scale_m[,2],col="blue")
  lines(type="b",scale_m[,3],col="black")
  lines(type="b",scale_m[,4],col="green")
  legend("topright",c("CV MSE","AIC","BIC","Cp"),
         col=c("red","blue","black","green"),lty=1)
}

plot_err(Rep_backward100[-c(5,6)],100)
plot_err(Rep_backward1000[-c(5,6)],1000)


plot_model_sizes = function(metrics,n){
  m = lapply(metrics, function(X) {apply(X,1,which.min)})
  par(mfrow=c(2,2),oma=c(0,0,3,0))
  hist(m$MSE,main="MSE",xlim=c(0,21),breaks=0:21,xlab="p")
  hist(m$AIC,main="AIC",xlim=c(0,21),breaks=0:21,xlab="p")
  hist(m$BIC,main="BIC",xlim=c(0,21),breaks=0:21,xlab="p")
  hist(m$Cp ,main="Cp ",xlim=c(0,21),breaks=0:21,xlab="p")
  title(paste("Empirical distribution of best model sizes\nfor n = ",n),outer=TRUE)
  par(mfrow=c(1,1),oma=c(0,0,0,0))
}

plot_model_sizes(Rep_backward100[-c(5,6)],n=100)
plot_model_sizes(Rep_backward1000[-c(5,6)],n=1000)



# Proportions -------------------------------------------------------------

## # left out first 15
prop_left = function(metrics){
  optimum = lapply(metrics[-c(5,6)], function(X) {apply(X,1,which.min)})
  proportions_left = optimum
  for (i in 1:B){
    proportions_left$MSE[i] = 1-length(intersect(metrics$variables[i,1:(optimum$MSE[i]-1)],metrics$names[i,1:15]))/15
    proportions_left$AIC[i] = 1-length(intersect(metrics$variables[i,1:(optimum$AIC[i]-1)],metrics$names[i,1:15]))/15
    proportions_left$BIC[i] = 1-length(intersect(metrics$variables[i,1:(optimum$BIC[i]-1)],metrics$names[i,1:15]))/15
    proportions_left$Cp[i] =  1-length(intersect(metrics$variables[i,1:(optimum$Cp[i] -1)], metrics$names[i,1:15]))/15
  }
  proportions_left
}

proportions_left100 = prop_left(Rep_backward100)
proportions_left1000 = prop_left(Rep_backward1000)


## # kept last 5
prop_kept = function(metrics){
  optimum = lapply(metrics[-c(5,6)], function(X) {apply(X,1,which.min)})
  proportions_kept = optimum
  for (i in 1:B){
    proportions_kept$MSE[i] = length(intersect(metrics$variables[i,1:(optimum$MSE[i]-1)],metrics$names[i,16:20]))/5
    proportions_kept$AIC[i] = length(intersect(metrics$variables[i,1:(optimum$AIC[i]-1)],metrics$names[i,16:20]))/5
    proportions_kept$BIC[i] = length(intersect(metrics$variables[i,1:(optimum$BIC[i]-1)],metrics$names[i,16:20]))/5
    proportions_kept$Cp[i] =  length(intersect(metrics$variables[i,1:(optimum$Cp[i])-1 ],metrics$names[i,16:20]))/5
  }
  proportions_kept
}

proportions_kept100 = prop_kept(Rep_backward100)
proportions_kept1000 = prop_kept(Rep_backward1000)



# Plot --------------------------------------------------------------------

plot_kept = function(proportions_kept,n){
  par(mfrow=c(2,2),oma=c(0,0,3,0))
  hist(proportions_kept$MSE,main="MSE",xlim=c(0,1),xlab="ratio")
  hist(proportions_kept$AIC,main="AIC",xlim=c(0,1),xlab="ratio")
  hist(proportions_kept$BIC,main="BIC",xlim=c(0,1),xlab="ratio")
  hist(proportions_kept$Cp, main="Cp",xlim=c(0,1),xlab= "ratio")
  title(paste("Proportion of times each coefficient was put in\nfor n = ",n),outer=TRUE)
  par(mfrow=c(1,1),oma=c(0,0,0,0))
}

plot_kept(proportions_kept100,100)
plot_kept(proportions_kept1000,1000)


plot_left = function(proportions_left,n){
  par(mfrow=c(2,2),oma=c(0,0,3,0))
  hist(proportions_left$MSE,main="MSE",xlim=c(0,1),xlab="ratio")
  hist(proportions_left$AIC,main="AIC",xlim=c(0,1),xlab="ratio")
  hist(proportions_left$BIC,main="BIC",xlim=c(0,1),xlab="ratio")
  hist(proportions_left$Cp, main="Cp",xlim=c(0,1),xlab= "ratio")
  title(paste("Proportion of times each coefficient was left out\nfor n = ",n),outer=TRUE)
  par(mfrow=c(1,1),oma=c(0,0,0,0))
}

plot_left(proportions_left100,100)
plot_left(proportions_left1000,1000)

# Methods which overfit:
# Methods which underfit:






