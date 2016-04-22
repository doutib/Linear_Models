load('family.rda')
attach(family)

# ## 1

intercept = rep(1,length(weight))
X = cbind(intercept,height, bmi)
betahat = solve(crossprod(X,X), t(X) %*% weight)
residuals = weight - X %*% betahat

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

# bmi is increasing with both weight and height


# ## 3








