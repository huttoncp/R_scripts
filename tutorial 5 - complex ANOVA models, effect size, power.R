#session 5: inferential statistics part 2: complex ANOVA models, effect size, power####

#effect coding setup####
options(contrasts=c("contr.sum","contr.poly")) #sets deviation (or sum-to-zero) coding for definition of effects, use before running any anovas or regressions

#import data####
library(car) #contains the capital "A" Anova function with type III SS and more detailed diagnostics.
#the car package also contains some useful practice data

OBdata <- OBrienKaiser #import the O'Brien-Kaiser repeated measures dataset 

str(OBdata) #review the structure of the data

#examine the data in a box plot####
boxplot(fup.2 ~ treatment*gender, data=OBdata)

#____2+ factorial ANOVA with categorical predictors____####

names(OBdata) #show just the column/variable names

#step 1: specify the linear model structure
lm.1 <- lm(fup.2 ~ treatment*gender, data=OBdata)

#step 2: analyze the model using ANOVA with type III SS
Anova(lm.1, type=3) # type III SS and F Tests

#simple main effects analysis#### 
#this is the easy way. Alternative is to do pairwise t-tests using subsets of the data. Ask me if interested.

library(emmeans)
emmeans(lm.1, ~ treatment | gender, contr="pairwise") #to evaluate simple main effects of A for each level of B use "|" between factors. p-value adjustment options include none, holm, tukey (default), and bonferroni
emmeans(lm.1, ~ gender | treatment, contr="pairwise") #to evaluate simple main effects of A for each level of B use "|" between factors. p-value adjustment options include none, holm, tukey (default), and bonferroni

#alternative ANOVA method which also gives generalized eta-squared####
library(afex)
OBdata$ID <- c(1:16) #add an ID variable

aov.1 <- aov_ez(id = "ID", dv = "fup.2", data = OBdata, between = c("treatment", "gender"))
summary(aov.1) #prints the ANOVA results summary table

#pariwise post-hoc comparisons####

emmeans(aov.1, ~ treatment*gender, contr="pairwise", adjust = "tukey") 
#other adjustment options include bonferroni (bonf), Holm (holm), or "none" (if justified by sig. main effect or interaction)
#tukey's adjustment is recommended when doing all possible pairwise comparisons

emmeans(aov.1, ~ treatment, contr="pairwise", adjust = "holm") 

emmeans(aov.1, ~ gender, contr="pairwise", adjust = "none") 

#____ANOVA diagnostics____####
outlierTest(lm.1) # Bonferonni p-value for most extreme obs based on the linear model residuals
#Change "n.max = " argument if you want more than 10 values to be shown (if there are more than 10 outliers)

outlierTest(lm.1, cutoff = Inf) #add "cutoff = Inf" argument to show all values ranked in order of "extremeness"

shapiro.test(residuals(lm.1)) #Shapiro-wilks test for normality for each group based on the residuals of your linear model, reject normality assumption if p < .05, may want to reject at .1 due to test being under powered

leveneTest(lm.1) #Levene's test for homogeneity of variance in models with 2 or more factors.

bartlett.test(fup.2 ~ treatment, data=OBdata) #test for homogeneity of variance if you have a single predictor
#reject normality assumption if p < .05, or maybe .1 as test is usually under-powered. 
#If data is non-normal use other tests

ncvTest(lm.1) #Breush-Pagan test for unequal variance in a linear model
#For more on homogeneity of variance see also http://www.cookbook-r.com/Statistical_analysis/Homogeneity_of_variance/

#visualizing devations from normality using a Quantile-Quantile (QQ) plot####
#https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot
#compares the x-y distrubtion of the model residuals 
#to the theoretical distribution of residuals if they had been obtained from a normal distribtion.
#If the linear model residuals are normally distributed, 
#then each point should fall near the regression line
#values that fall outside of the 95% confidence envelope represent 
#extreme values (potential outliers) that may be skewing the sample distribution.
lm.1 <- lm(pre.1 ~ treatment*gender, data=OBdata)

qqPlot(lm.1, distribution="norm", main="QQ Plot", ylab="Studentized Residuals", 
       col="black", col.lines="green", envelope=.95) #graphical estimate of normality with 95% confidence envelope

#additional methods of detecting influential observations: Cook's distance and Hat values####

windows()
infIndexPlot(lm.1, cex.lab = 1.5, cex = 2) 
summary(influence.measures(lm.1)) #identify potential outliers based on cook's D and Hat values

windows()
plot(lm.1, which=4) #plot of cook's D to identify potential influencial observations
(cutoff <- 4/(summary(lm.1)$df[2]-1)) #cutoff value for Cook's D is 4/(N-k-1) or 4/(denominator total df - 1)
abline(h=cutoff, col="red", lty=2) #add the cut-off line.

#model diagnostics all at once :)####
library(gvlma)

gvlma(lm.1, alphalevel = 0.05) 

#kurtosis and skewness are relevant to the normality assumption. 
#If your data are skewed, consider transforming the DV (e.g. log transform)

#heterscedasticity measure evaluates the equal (i.e. homogeneous) variance assumption

#global stat measure evaluates linearity of the model

#link function measure evaluates whether your dependent measure is continuous (expected)
#or if you should be using logistic regression or another method that requires a different link function

#if the above assumptions are not met consider a different test for more accurate results, e.g. a permutation test

#____ANCOVA: a mixture of categorical and continous predictors____####
data <- mtcars
str(data)

data$cyl <- as.factor(data$cyl) #convert number of cylinders to a factor
levels(data$cyl)

#let's say we want to know if engine size (# of cylinders) predicts miles per gallon
#independently of the weight of the vehicle, e.g. we want to evaluate wt as a co-variate

lm.1 <- lm(mpg ~ wt*cyl, data=data) #ANCOVA includes a continuous predictor
Anova(lm.1, type=3) # type III SS and F Tests

lm.1 <- lm(mpg ~ wt + cyl, data=data)
Anova(lm.1, type=3) # type III SS and F Tests


#____within-subjects ANOVA____####

library(car)
OBdata <- OBrienKaiser

str(OBdata)

#1st step: add an ID variable if there isn't one (needed for the error component of the repeated measures ANOVA model)####
OBdata$ID <- c(1:16)

#2nd step: convert a wide-format data set to long-format####
library(reshape2)
OBdata.long <- melt(data = OBdata, 
                    id.vars=c("ID", "treatment", "gender"), #variables you want to remain unchanged
                    measure.vars = c("fup.1", "fup.2", "fup.3", "fup.4", "fup.5"), #columns you want to combine 
                    variable.name = c("trial"), #new key variable for the old measurement column names
                    value.name=c("score")) #new value variable for the measurements that used to be spread across multiple columns
str(OBdata.long)

#alternative, more intuitive way using sequential functions from the dplyr and tidyr packages from the "tidyverse"
library(dplyr)
library(tidyr)

OBdata.long <- OBdata %>% #using OBdata
  select(ID, treatment, gender, fup.1:fup.5) %>% #select all the columns of interest. Can also use contains("text") or -contains("text") to include or exclude variables with "text" in the name
  gather(key = "trial", value = "score", #specify the names of the new key and value variables
                  c(fup.1:fup.5)) #Measurement columns to be gathered. Could instead use -c(ID, treatment, gender) to specify all other columns

str(OBdata.long)

OBdata.long <- OBdata %>% #using OBdata,
  select(ID, treatment, gender, contains("fup")) %>% #select all the columns of interest.
  gather(key = "trial", value = "score", #specify the names of the new key and value variables
         -c(ID, treatment, gender)) #Measurement columns to be gathered. 
str(OBdata.long)

#convert trial variable to a factor instead of character variable
OBdata.long$trial <- as.factor(OBdata.long$trial)

#3rd step: model specification and analysis####
#afex also includes effect size (ges) and diagnostic measure (sphericity assumption)

library(afex)
names(data)
?aov_car

rm.aov.1 <- aov_car(score ~ trial + Error(ID/trial), data = OBdata.long) #specify model, error term is key

rm.aov.1 #anova output with generalized eta squared effect estimates
#small effect ges = 0.01, medium effect ges = 0.06, and large effect ges = 0.14 (Olejnik & Algina, 2003)

summary(rm.aov.1) #Anova output with sphericity testing and GG/HF corrections, use corrected p-values if sphericity is violated

#____repeated-measures mixed (aka split-plot) ANOVA____####
library(afex)

names(OBdata.long)
rm.aov.1 <- aov_car(score ~ treatment*gender*trial + Error(ID/trial), data = OBdata.long) #specify model

rm.aov.1
#small effect ges = 0.01, medium effect ges = 0.06, and large effect ges = 0.14 (Olejnik & Algina, 2003)

#if you want to round values to 3 significant figures
round(rm.aov.1$anova_table, digits = 3) #anova output with generalized eta squared effect estimates

summary(rm.aov.1) #Anova output with sphericity testing info and GG/HF corrections, use corrected p-values if sphericity is violated

#post-hoc tests####
#tukey corrected pairwise comparisons for each level of the within subjects factor
#https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html#pairwise
#https://cran.r-project.org/web/packages/emmeans/vignettes/FAQs.html

library(emmeans)
emmeans(rm.aov.1, ~  treatment | trial, contr="pairwise", adjust = "tukey")

emmeans(rm.aov.1, ~  trial | treatment*gender, contr="pairwise", adjust = "holm")

levels(OBdata$treatment)

#polynomial contrasts####
(c.output <- emmeans(rm.aov.1, ~ trial | treatment*gender, contr="poly", adjust = "none"))


#pairwise comparisons between the polynomial contrasts for each group
contrast(c.output, method = "pairwise", by = "treatment") #pairwise comparisons of polynomial contrasts for each group with Tukey correction
contrast(c.output, method = "pairwise", by = "gender")
contrast(c.output, method = "pairwise", by = NULL) #pairwise comparisons of polynomial contrasts for each group with Tukey correction


#What if we want to know whether the linear change across trials differs between groups?
#i.e. does the linear change for males differ significantly from the linear change for females

#pairwise comparisons of the linear contrast for each group 
(output.df <- as.data.frame(c.output$contrasts)) #convert output to a data frame so that the row numbers are displayed

lin.output <- emmeans(rm.aov.1, ~ trial | treatment*gender, contr="poly")$contrasts[c(1,5,9,13,17,21),] #select just the linear output based on rows

#easier alternative method, use the grep function to find the row index numbers for you by matching a string pattern
(lin.rows <- grep(pattern = "linear", x = output.df$contrast))

(lin.output <- emmeans(rm.aov.1, ~ trial | treatment*gender, contr="poly")$contrasts[lin.rows,]) #select just the linear output based on rows

#pairwise comparisons split by each factor
contrast(lin.output, method = "pairwise", by = "treatment", adjust="tukey") #comparisons between genders, split by treatment group
contrast(lin.output, method = "pairwise", by = "gender", adjust="tukey") #comparisons between treatment groups, split by gender

#all comparisons
contrast(lin.output, "pairwise", by = NULL, adjust="tukey") #setting by = NULL removes the factor partitioning and enables interaction testing

#remove redundant text
all.output <- contrast(lin.output, "pairwise", by = NULL, adjust="tukey") #setting by = NULL removes the factor partitioning and enables interaction testing
all.output <- as.data.frame(all.output)
all.output$contrast <- gsub(pattern = "linear,", replace = "", x = all.output$contrast) #replace all instances of "linear," with blank text
all.output

#____MANOVA____####
#multiple dependent variables (not repeated measures)####

#for more info: http://online.sfsu.edu/efc/classes/biol710/manova/MANOVAnewest.pdf
#https://en.wikipedia.org/wiki/Multivariate_analysis_of_variance
#https://rpubs.com/aaronsc32/manova-test-statistics

#multivariate statistics####
library(car)
mlm.1 <- lm(cbind(fup.3, fup.4, fup.5) ~ treatment*gender, data = OBdata) #specify all DVs on left of ~ using the cbind function and all IVs on right of ~
Anova(mlm.1, type="3") #prints the pillai multivariate results, tells you if there is a difference between groups on one or more of the dependent measures

summary(Anova(mlm.1, type="3"), multivariate=T) #additional details and alternative multivariate tests, e.g. Wilks

#univariate follow-up tests####
summary(Anova(mlm.1, type="3"), multivariate=F, univariate=T) #get just the univariate test output

#multi- and uni-variate test results together####
summary(Anova(mlm.1, type="3"), multivariate=T, univariate=T) #multivariate and univariate results

#condensed combined output in 2 lines####
Anova(mlm.1, type="3") #prints the pillai multivariate results, tells you if there is a difference between groups on one or more of the dependent measures
summary(Anova(mlm.1, type="3"), multivariate=F, univariate=T) #get just the univariate test output

#____effect sizes for ANOVA____####

#eta and partial eta-squared measures of effect size for ANOVA using lsr package
#like R-squared for regression, eta-squarted describes the proportion of variability 
#associated with an effect when the variability associated with all other effects identified in the analysis 
#has been removed from consideration. 

#small effect = 0.01, medium effect = 0.06, and large effect = 0.14 for both estimates: Richardson, J. T. (2011). Eta squared and partial eta squared as measures of effect size in educational research. Educational Research Review, 6(2), 135-147.
#ref: Fritz, C. O., Morris, P. E., & Richler, J. J. (2012). Effect size estimates: current use, calculations, and interpretation. Journal of experimental psychology: General, 141(1), 2.
#tends to overestimate effects seen in small samples, therefore only use for larger samples (if you have a small sample use omega-squared instead)

#effect size for one-way or 2+ factorial ANOVA#### 
library(lsr)

lm.1 <- lm(fup.5 ~ treatment * gender, data=OBdata)

eta.sq <- etaSquared(lm.1, type=3) #eta-squared and partial eta-squared estimated using type III SS

etaSquared(lm.1, type=3, anova=T) #eta-squared and partial eta-squared with ANOVA table

#can also quickly get partial eta-squared manually via SSeffect/(SSeffect + SSresidual) or F(df1)/(F(df1)+df2)

Anova(lm.1, type=3) #confirm that the etaSquared fxn gives you an accurate anova table

#effect size for split-plot or mixed factorial ANOVA####
library(afex)
rm.aov.1 <- aov_car(score ~ treatment*gender*trial + Error(ID/trial), 
                    data = OBdata.long) #specify model

rm.aov.1 #anova output with generalized eta squared effect estimates,
#generalized eta squared is more appropriate than partial eta-squared for repeated-measures data
#small effect ges = 0.01, medium effect ges = 0.06, and large effect ges = 0.14 (Olejnik & Algina, 2003)

rm.aov.1 <- aov_car(score ~ treatment*gender*trial + Error(ID/trial), 
                    data = OBdata.long, 
                    anova_table = list(es = "pes")) #specify "pes" for partial eta squared as the effect size measure
rm.aov.1 #anova output with parial eta-squared effect estimates instead if you want them

rm.aov.1 <- aov_car(score ~ treatment*gender*trial + Error(ID/trial), 
                    data = OBdata.long, 
                    anova_table = list(es = c("ges", "pes"))) #can get both ges and pes!
rm.aov.1 #anova output with parial eta-squared and generalized eta-squared effect estimates instead if you want them

#____Association strength____#### 
# Olejnik S and Algina J (2003). Generalized Eta and Omega Squared Statistics: Measures of effect size for some common research designs. Psychological Methods 8(4) 434-447. 
# Lakens, D. (2013). Calculating and reporting effect sizes to facilitate cumulative science: a practical primer for t-tests and ANOVAs. Frontiers in psychology, 4, 863.

#use omega-squared and partial omega squared instead of eta-squared for small samples
library(sjstats)
lm.1 <- lm(fup.5 ~ treatment, data=OBdata)
omega_sq(lm.1, partial = F)
omega_sq(lm.1, partial = T)


#____power analysis____####
#requires cohen's f or f^2 for some functions
#very useful for planning experiment sample sizes using pilot data

#cohen's f (effect size)
lm.1 <- lm(fup.5 ~ treatment, data=OBdata)
summary(lm.1)
regr.1 <- summary(lm.1) #store regression summary object for adj.r.squared extraction below
sqrt(regr.1$r.squared/(1-regr.1$r.squared)) #cohen's f (a measure of effect size) for the overall model
#0.1 is a small effect, 0.25 is a medium effect, and 0.40 is a large effect.

#easy way using cohens_f() function from sjstats package
library(sjstats)
cohens_f(lm.1)

#correlation power####
library(pwr)
?pwr.r.test

#from the correlation and regression tutorial (#3)
library(dplyr)
data.mat <- as.matrix(select(mtcars, c(mpg, hp)))

library(Hmisc)
(cor.mat <- rcorr(data.mat, type=c("pearson"))) #pearson correlations and p-values for a specified matrix (x), 
#gives r-values, then # of observations (n) used, then p-values

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

#can also get the pearson r value (magnitude) from a regression table
lm.1 <- lm(mpg ~ hp, data=mtcars)
(regr.1 <- summary(lm.1))
(corr.value <- sqrt(regr.1$r.squared))


pwr.r.test(n = NULL, r = corr.value, sig.level = 0.05, power = 0.8) #fill in 3, get the last 1

#t-test power (equal variance)####
Mdata <- subset(OBdata, gender=="M")
Fdata <- subset(OBdata, gender=="F")

(sample.n <- length(Mdata$fup.5))

(dm <- (mean(Mdata$fup.5)-mean(Fdata$fup.5))) #get the difference between means for x1 and x2

(s <- sqrt((var(Mdata$fup.5) + var(Fdata$fup.5))/2)) #pooled standard deviation

power.t.test(n=sample.n, delta=dm, sd=s, sig.level=0.05, power = NULL, type="two.sample",
             alternative="two.sided") #evaluate the power of a t-test based on the sample size, sd, etc. could specify power instead of n to get the n required for a certain power level. 

#power calculation for welch t-test (unequal variance)####
y1 <- Mdata$fup.5
y2 <- Fdata$fup.5

(dm <- (mean(y1)-mean(y2))) #get the difference between means for x1 and x2

(sp.un <- sqrt(((length(y1)-1)*var(y1))+
                ((length(y2)-1)*var(y2))/
                (length(y1)+length(y2)-2))) #pooled standard deviation if the sample sizes are unequal

power.t.test(n = sample.n, d = dm, sd=sp.un, power= NULL, sig.level = 0.05, type = c("two.sample")) #specify the smaller n if sample sizes are unequal

#F-test (ANOVA) power####

#get cohen's f effect size metric 
lm.1 <- lm(fup.5 ~ treatment, data=OBdata)

library(sjstats)
(lm.1.f <- cohens_f(lm.1))
(treatment.f <- lm.1.f$cohens.f)

library(pwr)
pwr.anova.test(k=3, n = 8, f = treatment.f, sig.level=.05, power= NULL) #gives the number of subjects required to get a power of .8 if you have 3 groups, given a cohen's f of 0.28, and an alpha level of .05

pwr.anova.test(k=3, n = NULL, f = treatment.f, sig.level=.05, power= 0.8) #gives the number of subjects required to get a power of .8 if you have 3 groups, given a cohen's f of 0.28, and an alpha level of .05


#linear model power####
?pwr.f2.test

lm.1 <- lm(fup.5 ~ treatment, data=OBdata)
(regr.1 <- summary(lm.1))
regr.1$df
df.n <- regr.1$df[1]-1
df.d <- regr.1$df[2]

(lm.1.f <- sqrt(regr.1$r.squared/(1-regr.1$r.squared))) #cohen's f (a measure of effect size) for the overall model
#0.1 is a small effect, 0.25 is a medium effect, and 0.40 is a large effect.


pwr.f2.test(u = df.n, v = df.d, f2 = (lm.1.f^2), 
            sig.level = 0.05, power = NULL) # u and v are the numerator and denominator degrees of freedom, f2 is cohen's f-squared. only one part should be specified as NULL

pwr.f2.test(u = df.n, v = NULL, f2 = (lm.1.f^2), 
            sig.level = 0.05, power = 0.8) # u and v are the numerator and denominator degrees of freedom, f2 is cohen's f-squared. only one part should be specified as NULL

#effect size and power at the same time!####
library(sjstats) #contains the cohens_f function
library(car)
lm.1 <- lm(fup.4 ~ treatment*gender, data=OBdata)
anova_stats(Anova(lm.1, type = 3)) #cohen's f, power, and other useful values for each predictor

#purrr::map linear modelling####
library(tidyverse)
#linear models
EPM_data %>% 
  select_if(is.numeric) %>%
  map(~summary(lm(. ~ EPM_data$group)))