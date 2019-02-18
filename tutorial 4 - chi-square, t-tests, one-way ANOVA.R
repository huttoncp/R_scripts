#session 4: inferential statistics part 1: chi-square tests, t-test tests, one-way ANOVA####
#effect coding setup####
options(contrasts=c("contr.sum","contr.poly")) #sets deviation (or sum-to-zero) coding for definition of effects, use before running any anovas or regressions

#___CATEGORICAL VARIABLES####
#chi-square test####
#useful for comparing frequency/counts of values between two+ categorical variables
?chisq.test
#http://rcompanion.org/rcompanion/b_05.html
#https://www.r-bloggers.com/chi-squared-test/

#say we want to know if the incidence/frequency of self-reported memory impairment 
#1 month after a concussion (cause unknown) differs between older adults (>75) who exercise vs those who don't
#we survey 20000 individuals who exercise regularly and 30000 who do not (sedentary).
#of the exercisers, 18.6% report some degree of memory impairment
#while 29.8% of sedentary seniors report some degree of memory impairment
#we can compare these counts using a chi-square test (typically used for large samples)

#1. Entering the frequency/count data into vectors
incidence_exe <- round(20000*0.186, digits = 0) #compute the provided count, rounded to a whole number, for exercisers with memory problems.
incidence_sed <- round(30000*0.298, digits = 0) #digits argument determines number of decimal places to be used

exercise <- c(20000-incidence_exe, incidence_exe) #individuals without memory problems, individuals with memory problems
sedentary <- c(30000-incidence_sed, incidence_sed) #individuals without memory problems, individuals with memory problems

#2. combining the row vectors in matrix, then converting the matrix into a data frame (helpful for naming)
df <- as.data.frame(rbind(exercise, sedentary))

#3. assigning column names to this data frame
names(df) = c("normal", "impaired") #rename the columns

#4.display data and chi-square test
df
chisq.test(df, correct = F) #correct = F argument removes the default Yates' continuity correction for small sample sizes

#if you instead just have vectors of category labels
df.2 <- read.csv("label_vectors.csv")

head(df.2, n = 10)

#4.display data and chi-square test
chisq.test(x = df.2$treatment, y = df.2$improvement, correct = F) #correct = F argument removes the default Yates' continuity correction for small sample sizes

#save as a table (provides counts)
(df.2 <- table(df.2$treatment, df.2$improvement))
chisq.test(df.2, correct = F) #correct = F argument removes the default Yates' continuity correction for small sample sizes

#chi-square goodness of fit test: comparing an observed count to an expected count (1 sample)####
#say the indicence of self-reported memory impairments among adults >75 years of age is 20%

incidence_exe <- round(20000*0.186, digits = 0) #compute the provided count, rounded to a whole number, for exercisers with memory problems.
exercise <- c(20000-incidence_exe, incidence_exe) #individuals without memory problems, individuals with memory problems

chisq.test(x = exercise, #observed frequencies/counts 
           p = c(0.8, 0.2)) #expected proportions

#chi-square goodness of fit test for more than 2 groups####
#Create lists of values for this exercise
observed_distribution <- c(6, 7, 10)           #Number of observations in each group
expected_distribution <- c(1/3, 1/3, 1/3)       #Expected distribution across groups

#Run Chi Squared test
chisq.test(observed_distribution, p=expected_distribution)


#effect size for chi-square test####

#odds ratio = Pr(event) / Pr(not event) = P / 1-P###

#for odds ratio > 1: 1.68 = small effect , 3.47 = medium effect, and 6.71 = large effect
#for odds ratio < 1: 0.595 = small effect, 0.288 = medium effect, and 0.149 = large effect
#Chen, H., Cohen, P., & Chen, S. (2010). How big is a big odds ratio? Interpreting the magnitudes of odds ratios in epidemiological studies. Communications in Statistics—Simulation and Computation®, 39(4), 860-864.

library(epitools)
df
df.mat <- as.matrix(df)
oddsratio.wald(df.mat, conf.level = 0.95) #provides odds ratios for chi-square tests with 95% (or other) confidence intervals
#if you have a small sample, add "correction = T" argument

#can re-order rows to get inverse comparison
df.mat2 <- df.mat[c(2, 1),] #move row 2 to the 1st position and row 1 to the 2nd position
oddsratio.wald(df.mat2)

#can also just calculate inverse odds
1/oddsratio.wald(df.mat)$measure[2,] #inverse odds, in this case = odds of memory problems for exercisers vs. sedentary individuals


#graph the relative odds####
library(tidyr)
library(ggplot2)

df

#combine variables using gather function from tidyr package
df$life_style <- as.factor(row.names(df)) #add lifestyle grouping variable for graph using row names

df.long <- gather(data = df, key = memory_impairment, value = Y, -life_style) #
df.long

ggplot(data = df.long, aes(x = life_style, y = Y, fill = memory_impairment)) +
  geom_bar(stat = "identity", position=position_dodge(), colour="black") + 
  ylab("count") + xlab("life style") + 
  theme_bw(base_size = 16)

#___CONTINUOUS VARIABLES########
?colnames
?rownames
#import data####
data <- mtcars #store the mtcars data in a generic dataframe called "data" for convenience

#box plots####
names(data)

boxplot(hp ~ cyl, data=data, ylab="horse power", xlab="#cylinders", 
        main="hp by cyl boxplot") #generates a nonparametric box plot using dv and iv from data set "df". 
#Solid line is median, box edges are 25th and 75th percentiles

#t-tests####
#one-sample t-test####
#confidence intervals are included
names(data)
mean(data$mpg)

?t.test
t.test(data$mpg, alternative="two.sided", mu=100) 
t.test(data$mpg, alternative="less", mu=21)
t.test(data$mpg, alternative="greater", mu=100)

#many distributions available e.g. binomial, chisquare, F, t. 
#Can use them to find critical values for t-tests, F-tests, etc.
abs(qt(0.05/2, 40)) #bidirectional critical T value with a 2-tailed alpha of 0.05 and 40 df. abs() gives an absolute value with the sign removed.

#cohen's d effect size
s <- sd(data$mpg)
m <- mean(data$mpg)
(d <- ((m-100)/s)) #cohen's d = [sample mean - population mean (or chance value, enter manually)]/ standard deviation

#alternatively, d also = t/sqrt(length(x))
t.result <- t.test(data$mpg, alternative="two.sided", mu=100) 
d <- t.result$statistic/sqrt(t.result$parameter + 1)
names(d) <- c("d")
d #cohen's d

#two-sample t-test####
names(data)
levels(data$cyl)
data$cyl <- as.factor(data$cyl)

t.test(subset(data, cyl=="4")$hp, subset(data, cyl=="6")$hp, var.equal=TRUE) #t-test comparing the equivalence of groups x1 and x2 assuming equal variances, can specify alternative hypotheses using alternative = c("two.sided", "less", "greater")

t.test(subset(data, cyl=="4")$hp, subset(data, cyl=="6")$hp, var.equal=FALSE) #t-test comparing the equivalence of groups x1 and x2 assuming equal variances, can specify alternative hypotheses using alternative = c("two.sided", "less", "greater")

data4 <- subset(data, cyl=="4")
data6 <- subset(data, cyl=="6")
t.test(data4$hp, data6$hp, var.equal=FALSE) #t-test comparing the equivalence of groups x1 and x2 assuming equal variances, can specify alternative hypotheses using alternative = c("two.sided", "less", "greater")

#effect size for t-tests: cohen's d####

dm <- (mean(subset(data, cyl=="4")$hp)-mean(subset(data, cyl=="6")$hp)) #get the difference between means for x1 and x2

sp.en <- sqrt((var(subset(data, cyl=="4")$hp)+var(subset(data, cyl=="6")$hp))/2) #pooled standard deviation if both sample sizes are equal
sp.un <- sqrt(((length(subset(data, cyl=="4")$hp)-1)*var(subset(data, cyl=="4")$hp))+
                ((length(subset(data, cyl=="4")$hp)-1)*var(subset(data, cyl=="4")$hp))/
                (length(subset(data, cyl=="4")$hp)+length(subset(data, cyl=="4")$hp)-2)) #pooled standard deviation if the sample sizes are unequal

d <- dm/sp.un # cohen's d. small effect = 0.2, medium effect = 0.5, large effect = 0.8 (Cohen, 1969)
d

#____MIX of CAT and CONT variables##### 

#one-way anova####
library(car)
OBdata <- OBrienKaiser
names(OBdata)
str(OBdata)
levels(OBdata$treatment)

#using basic R distribution of packages (installed with R)####
aov.1 <- aov(fup.2 ~ treatment, data = OBdata) #uses type I SS, only valid for a single categorical predictor
summary(aov.1)

1-pf(6.136, df1=2, df2=13) #calculate the p-value associated with an observed F statistic.

#alternative method which also gives generalized eta-squared####
library(afex)
OBdata$id <- c(1:16) #add an ID variable if needed
names(OBdata)
?aov_ez
aov.1 <- aov_ez(id = "id", dv = "fup.5", data=OBdata, between = c("treatment"))
aov.1
aov.1$anova_table


#post-hoc testing (if the factor has 3 or more levels)####
#do t-tests for a small number of comparisons or to break down interactions, or ANOVA on subsets of data
#install.packages("agricolae")
library(agricolae)

LSD.test(aov.1, "treatment", p.adj="none", group = F, console=T) #only use if you are doing 3 comparisons, if doing more, use Tukey's HSD or use p.adj="bonferroni"

TukeyHSD(aov.1, which="treatment") #all pairwise comparisons with Tukey corrected p-values

TukeyHSD(aov.1, which="treatment", conf.level=0.90) #Tukey's HSD test for pairwise comparisons at each level of B with a 90% confidence interval

pairwise.t.test(OBdata$fup.3, OBdata$treatment, p.adj = "bonf") #multiple pairwise planned comparisons with bonferroni adjustment where x is the response vector and g is the grouping vector
pairwise.t.test(OBdata$fup.3, OBdata$treatment, p.adj = "holm") #holm's sequential bonferroni test for planned pairwise comparisons

#alternative method using emmeans (generalizes to more complex models and contrasts)
library(emmeans)
emmeans(aov.1, ~ treatment, contr="pairwise", adjust="tukey") #to evaluate simple main effects of A for each level of B use "|" between factors. p-value adjustment options include none (lsd), holm, tukey (default), and bonferroni


#correcting p-values for multiple comparisons####

alpha.pc <- .05 #code per comparison alpha level
C <- 2 #code the number of comparisons you did, in this case "2"
(alpha.fw <- (1 - (1-alpha.pc)^C)) #determine the family wise alpha level
(alpha.pc <- alpha.fw/C ) #get the pc alpha to use for a given number of comparisons while holding fw alpha at .05 bonferroni style

p <- c(0.03,0.04) #create variable p composed of the p-values you are want an adjustment on

L <- length(p) #the number of comparisons you are doing

p.adjust(p, "holm", L) #does the holm p-value adjustment (less conservative)

p.adjust(p, "bonferroni", L) #does the bonferonni p-value adjustment (more conservative)

?p.adjust


