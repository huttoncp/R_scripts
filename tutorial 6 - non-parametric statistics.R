#session 6: non-parametric statistics####

#some options if there are outliers or your data are not normally distributed
#and/or are heteroskedastic (unequal variance between groups)
#although less powerful than parametric tests if the assumptions of those tests are met,
#these methods are more powerful (and thus preferable) than parametric tests if the assumptions are not met.
#see the tutorial 4 R script for diagnostic tests.
#there are some other tests available for R but not covered here due to time constraints, e.g. mann-whitney U test
#see https://www.statmethods.net/stats/nonparametric.html if you are interested.

#1. load or import your data
data <- mtcars #store the mtcars data in a generic dataframe called "data" for convenience

names(data)
levels(data$cyl)
class(data$cyl)

data$cyl <- as.factor(data$cyl)

OBdata <- OBrienKaiser

#2. take a look at it
boxplot(hp ~ cyl, data = data) #boxplots are a useful 1st step in seeing if there are outliers or the data may be skewed

#3. fit a model and run diagnostic tests -> see tutorial 4 script.

####################################################################################################
#welch t-test for samples with unequal variance####
#if data are also non-normally distributed use a permutation test instead

levels(data$cyl)
data$cyl <- as.factor(data$cyl)

t.test(subset(data, cyl=="4")$hp, subset(data, cyl=="6")$hp, var.equal=TRUE) #t-test comparing the equivalence of groups x1 and x2 assuming equal variances, can specify alternative hypotheses using alternative = c("two.sided", "less", "greater")

t.test(subset(data, cyl=="4")$hp, subset(data, cyl=="6")$hp, var.equal=FALSE) #t-test comparing the equivalence of groups x1 and x2 assuming equal variances, can specify alternative hypotheses using alternative = c("two.sided", "less", "greater")

#purrr::map'd t-tests####
cyl4.6_data <- data %>% filter(cyl =="4" | cyl =="6") #can only compare 2 sample groups with t-test function

library(tidyverse) #contains dplyr and purrr packages, could load separately using library(dplyr) & library(purrr)

cyl4.6_data %>% 
  select_if(is.numeric) %>%
  map(~t.test(. ~ cyl4.6_data$cyl, var.equal=T)) #substitute "." for y, since map will loop over each "y" column


#Welch equivalent to a one-way (single-factor) ANOVA####
oneway.test(hp ~ cyl, data = data) #Welch correction for a one way anova assuming no homogeneity of variance


#Kruskal-Wallis test alternative to one-way ANOVA####
#for simple models where normality and homogeneity of variance are both violated
kruskal.test(data$hp, data$cyl) # kruskal-wallis test for variance assuming no homogeneity of variance or normality
?kruskal.test

#2-way ANOVA using trimmed means (ignores extreme points AKA outliers before calculating the means)####
library(WRS2)
?t2way
t2way(fup.5 ~ treatment*gender, data = OBdata) #add tr = 0.1 or other number to change the trim level for the mean
mcp2atm(fup.5 ~ treatment*gender, data = OBdata) #post-hoc comparisons using trimmed means

#2-way ANOVA using medians####
?med2way
med2way(fup.5 ~ treatment*gender, data = OBdata) 

#post-hoc tests
mcp2a(fup.5 ~ treatment*gender, data = OBdata, est = "median") #post hoc tests

####################################################################################################
#K-S tests to compare the group distributions non-parametritcally (need to correct p-values for multiple comparisons)####
library(stats)
?ks.test() #lowercase
names(OBdata)
levels(OBdata$gender)
Mdata <- subset(OBdata, gender=="M") #extracts the subset of the dataframe for group A of variable x and stores it in a new data frame called xA
Fdata <- subset(OBdata, gender=="F") #extracts the subset of the dataframe for group B of variable x and stores it in a new data frame called xB

ks.test(Mdata$post.3, Fdata$post.3) # can specify alternative hypotheses as usual, ex. alternative="greater", then tests for xA > xB
####################################################################################################


####################################################################################################
#robust regression and ANOVA#####
#good if outliers are a problem. Also good if heterskedasticity is an issue (unequal variance)

#using the rlm function which employs iterative re-weighted least squares using robust MM-estimators
#preferred for heteroskedastic data or data with outliers.
#assumes asymptotic normal distribution in majority of data aside from a few outliers. 
#If data are far from normally distributed or you have a very small sample use a resampling method (permutation test or bootstrapping) instead 
#using robust package can also get a measure of model bias/reliability

library(MASS)
library(car)
?rlm
rlm.1 <- rlm(fup.5 ~ treatment*gender, method="MM", data=OBdata)
summary(rlm.1) #regression table
Anova(rlm.1, type=3) #anova table from car package with p-values
?rlm


library(robust)
??robust
?lmRob
lm.1 <- lm(fup.5 ~ treatment*gender, data=OBdata)
summary(lm.1) #ordinary least-squares regression output
summary(lmRob(lm.1)) #robust regression output with test for bias

####################################################################################################

#resampling methods: permutation tests and bootstrapping####

#NOTES
#bootstrapping is primarily used to estimate the precision of your statistic via confidence intervals
#permutation tests are primarily used for hypothesis testing (obtaining p-values)

#both assume that your data represent a random sample and are therefore likely representative of the population
#both are most useful when you have small samples (n < or = 12)
#for experimental treatments, requires that subjects be randomly assigned to treatment conditions (usually the case anyways)
#by "resampling"(bootstrapping) or randomly "re-assigning" subjects to groups many times,
#you can get a more accurate estimate of the
#due to the resampling process, there is no requirement for the data to conform to a specific (e.g. normal) disribution
#the resampling or reassignment process is also robust to outliers.
#Note: a disadvantage of these methods is that they are very computationally intense for larger samples 
#and can take a while to run on a slower computer (less of an issue now). 
#Howewer, a newer PC or laptop can do 10,000,000 permutations in < 10 seconds. A recommended minimum is ~100,000.
#This was a serious issue when they were first developed which is why they are not used as frequently

####################################################################################################

#permutation tests#### 
#non-parametric tests for small samples, abnormally distributed, you have outliers, or unbalanced data (unequal group sizes)
#use a permutation test instead of a t-test or ANOVA if your group sample size is less than 6

library(coin)
oneway_test(fup.5 ~ treatment, data=OBdata, distribution=approximate(B=100000)) #permutation test alternative to single factor anova tha approximates an exact p-value based on B iterations, useful for small samples

library(perm)
permKS(fup.5 ~ treatment, data=OBdata, exact=T) #alternative to above method

permTS(Mdata$fup.5, Fdata$fup.5, alternative = "greater") # permutation alternative to a t-test to compare 2 groups with direction hypothesis options
?permTS

#purrr::map permutation tests####
names(EPM_data)
unique(EPM_data$Animal.ID)

EPM_data %>% 
  select_if(is.numeric) %>%
  map(~permTS(. ~ EPM_data$group))



####################################################################################################

#for complex models use lmPerm package#### 
#lmPerm has to be installed from github using devtools package

install.packages("devtools")#required to install R packages from github instead of CRAN
library(devtools)
install_github("mtorchiano/lmPerm")
library(lmPerm)
?aovp

#use aovp function for anovas and lmp function for linear models 
library(lmPerm)
aovp.1 <- aovp(fup.5 ~ treatment*gender, seq=F, maxIter = 1000000, Ca=0.0000001, data = OBdata) #anova model for permutation test using type III SS
summary(aovp.1) #prints the anova table with estimated p-values based on the random sample of possible permutations

#permutation anova effect size: partial eta-squared using formula SSeffect / (SSeffect+SSresidual)
summary.1 <- summary(aovp.1)
SS <- cbind(summary.1[[1]]$`R Sum Sq`)
p.eta.sq <- format((p.eta.sq <- SS/(SS+(SS[length(SS),]))), scientific=F)#partial eta-squared effect sizes for aovp permutation effects using formula SSeffect/SSresidual
colnames(p.eta.sq) <- c("partial eta-squared")
(perm.effects <- cbind(rownames(summary.1[[1]]), p.eta.sq))[-8,] #effect size table

library(lmPerm)
lmp.1 <- lmp(fup.5 ~ treatment*gender, seq=F, maxIter = 1000000, Ca=0.0000001,  data=OBdata) #linear model for permutation test using type III SS
summary(lmp.1) #regression output

#use aovp function for repeated-measures anova
OBdata <- OBrienKaiser
names(OBdata)

OBdata$ID <- c(1:16) #add a subject ID variable if there isn't one

library(reshape2)
OBdata.long <- melt(OBdata, id.vars=c("ID", "treatment", "gender"), 
                    measure.vars = c("pre.5", "post.5", "fup.5"), variable.name = c("trial"), 
                    value.name=c("score"))

aovp.1 <- aovp(score ~ treatment*gender*trial + Error(ID/trial), 
               seq=F, maxIter = 1000000, Ca=0.0000001, data = OBdata.long) #anova model for permutation test using type III SS
summary(aovp.1) #prints the anova table with estimated p-values based on the random sample of possible permutations

##########################################################################################################
#bootstrapping####

#using the simpleboot package####
library(simpleboot)

#one sample
median(mtcars$hp) # the sample median
boot.1 <- one.boot(mtcars$hp, median, R = 10000) #data, function name, R = #re-samples, add student = T if you want a studentized bootstrap
boot.ci(boot.1) #BCa gives percentile estimates for the CI but they are also adjusted for bias and skewness in the data
hist(boot.1) #red line is the observed value of the statistic

?one.boot

#difference between a statistic for two samples, e.g. median difference
Mdata <- subset(OBdata, gender=="M")
Fdata <- subset(OBdata, gender=="F")
boot.2 <- two.boot(Mdata$pre.1, Fdata$pre.1, median, R = 10000, student = F) #sample1$y, sample2$y, function name, R = # of resamples 
boot.ci(boot.2)
hist(boot.2)
median(Mdata$pre.1) - median(Fdata$pre.1)

#bootstrapping a correlation coefficient
library(simpleboot)
cor.boot <- pairs.boot(data$y1, data$y2, cor, R = 10000, student = F)
boot.ci(cor.boot)
hist(cor.boot)
?boot.ci

#bootstrapping a linear model: regression coefficient
lm.1 <- lm(mpg ~ wt, data = mtcars)
lboot <- lm.boot(lm.1, R=10000)
summary(lboot)

#vs. regular regression output
summary(lm.1) #ordinary least-squares regression output

#using the boot package
library(boot)

data <- mtcars

my.median = function(x, indices) {
  return(median( x[indices] ) )
}

hp.boot <-boot(data$hp, my.median, 10000)
hp.boot
boot.ci(hp.boot)
hist(hp.boot)
?hist

#R also has other options for bootstrapping, see https://www.statmethods.net/advstats/bootstrapping.html
