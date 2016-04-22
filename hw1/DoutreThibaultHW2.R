load('family.rda')
attach(family)

# ## 1

intercept = rep(1,length(weight))
X = cbind(intercept,height, bmi)
betahat = solve(crossprod(X,X), t(X) %*% weight)
betahat
residuals = weight - X %*% betahat
residuals

# ## 2

regcoef = function(df=family[4:5]){
  x = df[,1]
  y = df[,2]
  b = sum((y-mean(y))*(x-mean(x)))/sum((x-mean(x))^2)
  a = mean(y) - b*mean(x)
  return(c(a,b))
} 
regline = function(df=family[4:5]){
  coeff = regcoef(df)
  plot(height,weight, main = "Least Squares Regression Line",cex = bmi*0.06)
  abline(coeff[1],coeff[2])
}
regline()


# ## 3

library(rgl)

par3d(params=list(
  windowRect=c(100,100,600,600)))
view3d( theta = 70, phi = 20)

plot3d(X[,2],X[,3],weight,"height","bmi","weight",col='darkgray')
planes3d(betahat[2],betahat[3],c=-1,betahat[1],add=TRUE,alpha=.5,col="lightblue")
rgl.postscript("plot.pdf","pdf")

# ## 4

M = X[,1:2]
N = X[,3]
Y = weight

# step 1
gamma_1 = solve(crossprod(M,M), crossprod(M,Y))
gamma_1
f = Y - M %*% gamma_1

# step 2
gamma_2 = solve(crossprod(M,M), crossprod(M,N))
gamma_2
g = N - M %*% gamma_2

# step 3
gamma_3 = sum(f*g)/sum(g*g)
gamma_3
e = f-g*gamma_3

# step 4
betahatBis = c(gamma_1-gamma_2*gamma_3,gamma_3)
betahatBis

# Check validity 
abs(betahat-betahatBis)
