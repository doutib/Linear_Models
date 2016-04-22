###  HW 1   Your Name Here
###
###  Part 1.
###
### 1. put your code here for loading the data, write comments if anything
### isn't obvious.

### 2. Write the function regcoef() as described in the assignment.
### You can change what I did but this is a suggestion as to how I'd
### start off writing the function.  For example it's not necessary to
### have a default for df but if you are going to be
### running and testing your function it will save typing.

regcoef=function(df=family[4:5]){
  # This function takes input of a data frame and returns regression 
  # coefficients predicting the second variable from the first.
}

### 3. Similarly, write the regline() function.

regline=function(df=family[4:5]){
  # This function plots points and regression line, doesn't need to return
  # anything.
  coefs=regcoef(df)
}