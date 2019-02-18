#session 3: correlation, scatterplots, regression, and effect sizes#####

#specification of default contrasts in model comparisons for inferential stats#####
options(contrasts=c("contr.sum","contr.poly")) #sets deviation (or sum-to-zero) coding for definition of effects, use before running any anovas or regressions

#R and most packages have built in datasets that can be used to test functions, e.g. the motor trend car road test dataset called "mtcars"

data <- mtcars #store the mtcars data in a generic dataframe called "data" for convenience

#....correlations...............................................................####
names(data)
?mtcars
str(data) #show the structure of the data

names(data) #show the column/variable names

data.mat <- as.matrix(mtcars) #convert the data to matrix form

#subset using column index numbers####
data.mat <-as.matrix(data[,1:4]) #extracts columns of the data that you want correlations for (excludes the string factors)
data.mat <-as.matrix(data[,c(1, 3, 6)]) #extracts just the columns of interest from all data and stores them in a matrix
data.mat

#subset using column names, the select() function####
library(dplyr)
?select()
data.mat <- as.matrix(select(data, mpg:hp)) #use "-" to exclude columns, e.g. -mpg
data.mat <- as.matrix(select(data, c(mpg, hp, disp)))

#correlations####
library(Hmisc)
cor.mat <- rcorr(data.mat, type=c("pearson")) #pearson correlations and p-values for a specified matrix (x), 
#gives r-values, then # of observations (n) used, then p-values
cor.mat

######################################################################
#custom function to convert rcorr matrix to column for easier reading# 
flattenCorrMatrix <- function(cormat, nmat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    n = nmat[ut],
    p = pmat[ut])
}
#######################################################################
names(cor.mat) # provides the list of variables(columns) in your rcorr output matrix

cor.mat.v2 <- flattenCorrMatrix(cor.mat$r, cor.mat$n, cor.mat$P)
cor.mat.v2

#....scatterplots.................####
library(car) #contains the scatterplot function and many others - e.g. can add marginal boxplots.

names(mtcars)

scatterplot(mpg ~ hp, data = mtcars, smooth = T) #scatterplot with regression line. 
#set smooth = T to add lowess smooth line 
#boxplots can be turned off also. See ?scatterplot help page for more options like colours
?scatterplot
windows() #opens a plotting window

names(data)
scatterplot(mpg ~ disp | cyl, data = data, boxplots=T, smooth = F) #split the scatterplot by a grouping factor

library(car)

#3D scatterplots####
scatter3d(mpg ~ hp + cyl, data = data) #3D scatterplot


# 3D Scatterplot with Coloring and Vertical Lines
# and Regression Plane 
library(scatterplot3d) 
?scatterplot3d

s3d <- scatterplot3d(x = data$wt, y = data$disp, z = data$mpg, #have to specify dataset$variable since the function lacks a data argument
                     xlab = "weight", ylab = "displacement", zlab = "miles per gallon",
                     pch=16, highlight.3d=TRUE,
                     type="h", main="3D Scatterplot")

#add a linear regression plane
fit <- lm(mpg ~ wt+disp, data = data) 

s3d$plane3d(fit)

#scatterplot matrices####
names(data)
data.sub <- data[,c(1, 2, 4, 6)]
str(data.sub)

pairs(data.sub) #provides a basic scatterplot matrix

windows()
scatterplotMatrix(~ mpg + hp + wt | cyl, data=data.sub, main="Three Cylinder Options") #car scatterplot matrix

#linear regression####
data <- mtcars
names(data)

lm.1 <- lm(mpg ~ hp, data=data) #specifies the DV and IV to be compared in the linear model

summary(lm.1) #regression table

dummy.coef(lm.1) #gives all of the regression coefficients for the predictor(s)

confint(lm.1, level = 0.95) #provides confidence intervals for a linear model

#Multiple linear regression####
lm.1 <- lm(formula = mpg ~ hp + cyl, data = data) # formula = Y ~ X1 + X2. 
lm.1 # see the basic coefficients of your regression

summary(lm.1) #regression summary table

#evaluating predictor interactions
lm.1 <- lm(formula = mpg ~ hp + cyl + hp:cyl, data = data) #formula = Y ~ X1 + X2 + X1:X2, where X1:X2 represents the interaction between X1 and X2

lm.1 <- lm(formula = mpg ~ hp*cyl, data = data) # formula = Y ~ X1*X2. (X1*X2) is equivalent to (X1 + X2 + X1:X2)

summary(lm.1) #regression summary table
names(data)

#....effect size..................................................................####

#R^2 = the proportion of the dependent vairable explained by the independent variables####

#full model adjusted R^2 and partial R^2 values for each predictor and interaction
library(rsq) #contains the rsq function for obtaining the full and partial R-squared values for your linear model

rsq(lm.1, adj = F) #full model un-adjusted R^2.

rsq(lm.1, adj = T) #full model adjusted R^2.

summary(lm.1) #can also get the full model R^2 from the summary(lm) function, AKA the regression summary table

rsq.partial(lm.1, adj = T) #obtain partial R-squared values for each predictor. Add adj=T for adjusted R^2 values or adj = F for un-adjusted R^2.

#cohen's f for the overall linear model (using formula)####
regr.1 <- summary(lm.1) #store regression summary object for adj.r.squared extraction below

sqrt(regr.1$adj.r.squared/(1-regr.1$adj.r.squared)) #cohen's f (a measure of effect size) for the overall model
#0.1 is a small effect, 0.25 is a medium effect, and 0.40 is a large effect.


#cohen's f for each predictor####
library(sjstats) #contains the cohens_f function

?cohens_f
cohens_f(lm.1) #cohen's f: 0.1 is a small effect, 0.25 is a medium effect, and 0.40 is a large effect).

#cohen's f is also needed for power calculations (see tutorial 4 script)

#association strength####

# omega^2 of .01 is a small/weak association, .06 is a medium/moderate association, and .14 is a large/strong association

#omega^2####
library(sjstats)
omega_sq(lm.1, partial = F)

#partial omega^2####
omega_sq(lm.1, partial = T)

