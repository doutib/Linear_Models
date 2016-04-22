###  HW 1   Doutre Thibault

load('family.rda')
attach(family)

# ## 1

# Least Squares Regression Line
# y = b*x + a
b = sum((weight-mean(weight))*(height-mean(height)))/sum((height-mean(height))^2)
a = mean(weight) - b*mean(height)

plot(height,weight, main = "Least Squares Regression Line")
abline(a,b)

# ## 2

regcoef = function(df=family[4:5]){
  x = df[,1]
  y = df[,2]
  b = sum((y-mean(y))*(x-mean(x)))/sum((x-mean(x))^2)
  a = mean(y) - b*mean(x)
  return(c(a,b))
} 

regcoef()

# ## 3

regline = function(df=family[4:5]){
  coeff = regcoef(df)
  plot(height,weight, main = "Least Squares Regression Line")
  abline(coeff[1],coeff[2])
}

regline()


