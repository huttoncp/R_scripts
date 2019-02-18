#session 2: generating random numbers, exporting data, manipulating/reorganizing data, descriptive statistics, confidence intervals##########
setwd("C:/Users/CPH/Desktop/documents/McMaster/stats/R training Nov 2018/")

################################################################################################################
#random samping####
? sample

#obtaining a random sample and exporting it into a csv file for easy copying into data sheets
x <- c(1:32) #the range of subject ID numbers
x
x1 <- as.matrix(sample(x, 32)) #obtain random subject order and store it as a matrix
write.csv(x1, file = "random sample.csv") #export to csv file for easy copying, file will appear in your working directory

#obtaining a random sample for all subjects except 1 (ex if a subject dropped from the study)and exporting it into a csv file for easy copying into data sheets
x <- 201:232 #the range of subject ID numbers
length(x)
x1 <-
  as.matrix(sample(x[x != 229], 31)) #obtain random subject order skipping number 229 and store as a matrix
write.csv(x1, file = "mouse order.csv") #export to csv file for easy copying, file will appear in your working directory

#obtaining random subject order####
x <- 201:232
d1 <-
  sample(x, 32) #obtain random subject order and store as a matrix
d2 <-
  sample(x, 32) #obtain random subject order and store as a matrix
d3 <-
  sample(x, 32) #obtain random subject order and store as a matrix
d4 <-
  sample(x, 32) #obtain random subject order and store as a matrix
d5 <-
  sample(x, 32) #obtain random subject order and store as a matrix
d6 <-
  sample(x, 32) #obtain random subject order and store as a matrix
d7 <-
  sample(x, 32) #obtain random subject order and store as a matrix

week <-
  cbind(d1, d2, d3, d4, d5, d6, d7) #combine into single matrix
colnames(week) <-
  c("day1", "day2", "day3", "day4", "day5", "day6", "day7") #label the columns
week #show the new matric on the screen
write.csv(week, file = "mouse order.csv") #export to csv file for easy copying
? write.csv

write.csv(week, file = "mouse order.csv", row.names = FALSE) #add "row.names = FALSE" to get rid of the row names/index #s

getwd()
#reshaping data using reshape2####
#converting from wide to long format using the melt function from the reshape2 package
library(reshape2)
library(car)

? OBrienKaiser
OBdata <- OBrienKaiser
str(OBdata)
names(OBdata)

#add an ID or subject number variable
OBdata$ID <- c(1:16)

#if you only have a few columns to merge
? melt
data_long <-
  melt(
    OBdata,
    id.vars = c("ID", "gender", "treatment"),
    # ID variables - all the variables to keep but not split apart on
    measure.vars = c("pre.1", "pre.2", "pre.3"),
    # The source columns you want to combine
    variable.name = "pretest.hour",
    value.name = "measurement",
    na.rm = TRUE
  ) # Name of the destination column that will identify the original # column that the measurement came from

str(data_long)
levels(data_long$pretest.hour) #list the groups of the factor of interest, as ordered by R

levels(data_long$pretest.hour) <-
  c("h1", "h2", "h3") #change the group names, order matters!

levels(data_long$pretest.hour)[levels(data_long$pretest.hour) == "h1"] <-
  "hour1" #change the name of a specific group, order doesn't matter
levels(data_long$pretest.hour)[levels(data_long$pretest.hour) == "pre.2"] <-
  "hour2"
levels(data_long$pretest.hour)[levels(data_long$pretest.hour) == "pre.3"] <-
  "hour3"

levels(data_long$pretest.hour) #list the groups of the factor of interest, as ordered by R

# if you have many columns to merge don't specify measure.vars
data_long <-
  melt(
    OBdata,
    id.vars = c("ID", "treatment", "gender"),
    # ID variables - all the variables to keep but not split apart on (ie. the source columns you want to combine)
    variable.name = "pretest.hour",
    value.name = "measurement"
  )
class(data_long)

#converting from long to wide format using the dcast function also from the reshape2 package
library(reshape2)
data_wide <-
  dcast(
    data_long,
    ID + treatment + gender ~ pretest.hour,
    fun.aggregate = mean,
    value.var = "measurement"
  )  #"subject" and "sex" are columns we want to keep the same. "condition" is the column that contains the names of the new column to put things in. "measurement" holds the measurements

# data_tidy <- dcast(data_proteins,
#     subject + sex + experimental_group ~ accession_number, #vars you want to leave unchanged ~ key variable you want to spread from
#     fun.aggregate = mean,
#     value.var = "normalized_intensity")  


data_wide <-
  dcast(data_long, id ~ pretest.hour, mean, value.var = "measurement")

? dcast
str(data_wide)

#convert back to long format####
data.v2 <-
  melt(
    wide.data,
    id.vars = c("subject", "group", "genotype", "sex", "diet"),
    measure.vars = c("t21", "t22", "t23", "t24"),
    variable.name = c("trial"),
    value.name = c("time")
  )

#merging data####
#see: https://www.statmethods.net/management/merging.html

d1d2 <- rbind(d1, d2) #combine data frames vertically (by row)

#recode/label values####
#using dplyr, convert multiple values and turn unspecified ones into NAs
library(dplyr)
recode(x, '1' = "1st", '2' = "2nd") #variable/column name/number
#recoding of hrs per day
x <- c(1, 1, 2, 2, 2, 5, 1, 4, 2, 6, 67, 8)

#using base R coding, will change only the values specified
x[x == "1"] <- "first"
x[x == "1"] <- "100000"

#rename a variable or column####
#using dplyr
library(dplyr)
rename(OR, odds_ratio = OR) #rename a column: rename(data.frame, new.name = old.name)

#using base R code
names(data_wide)[names(data_wide) == "cond1"] <- "first"
names(data_wide)[names(data_wide) == "cond2"] <- "second"


#renaming factor levels####
#using dplyr
# Use recode_factor() to create factors with levels ordered as they are specified
x <- as.data.frame(cbind(1:80))
x$v2 <- as.factor(letters[1:4])
levels(x$v2)
str(x)
x$v2 <- recode_factor(x$v2, `c` = "cookies", `b` = "beer")
levels(x$v2)

#using base R code
levels(x$v2)[levels(x$v2) == "beer"] <- "B"
levels(x$v2)

#reordering factor levels (useful for graphing)####
#using forcats
#https://blog.rstudio.com/2016/08/31/forcats-0-1-0/
#fct_relevel() is similar to stats::relevel() but allows you to move any number of levels to the front.

sizes <- factor(c("small", "large", "large", "small", "medium"))
levels(sizes)

library(forcats)
? fct_relevel()

sizes <-
  fct_relevel(sizes, "small", "medium", "large") #can specify the position to insert the level after
levels(sizes)

sizes <-
  fct_relevel(sizes, "small", after = 2) #can specify the position to insert the level after
levels(sizes)



#changing the order of columns####
names(data_wide)
data_wide <-
  data_wide[, c(1, 2, 3, 6, 5, 4)] # assuming you have 5 columns, the 5th column will now be in between the 2nd and 3rd columns


############################################################################################################
#Basic descriptive statistics and data summaries####
mean(x$V1, na.rm = T) #mean of variable 'x'
sd(x$V1) #standard deviation of variable 'x'
summary(OBdata) #lists some basic stats for the data frame 'df'

max(x$V1) #maximum
median(x$V1) #median
min(x$V1) #minimum
sum(x$V1) #sum
var(x$V1) #variance
mean(x, na.rm = TRUE) #excludes NA values from analyses, in this case calculation of the mean

names(OBdata)
mean(OBdata$pre.2) #calculates the mean of variable 'x' within data frame 'df'

? tapply
tapply(OBdata$pre.1, OBdata$gender, mean) #provides a calcualtion of "statistic"(mean, sd, length etc) for the dependent variable "dv" and the independent variable "iv"

with(OBdata, tapply(pre.1, gender, mean)) #gives the mean for the scores(dv) in each group(iv), using data from the dataframe "data"

library(psych) #contains the describeBy function
data <- mtcars
? mtcars
describe(data$mpg) #gives descriptive stats for a variable of interest
describeBy(OBdata$pre.1, OBdata$treatment, mat = TRUE) #gives descriptive stats for each group of a data frame of interest.

library(Rmisc) #contains the summarySE and summarySEwithin functions
summarySE(
  data = OBdata,
  measurevar = "pre.1",
  groupvars = c("treatment"),
  na.rm = T
) #means, std dev, SE and CI for a DV split into groups on an IV. To get the upper and lower bounds of the CI just add or subtract the CI value from the mean

#confidence intervals with plotting option#####
library(Publish) #contains the ci.mean function and CI plotting function
? ci.mean
(ci1 <-
    ci.mean(OBdata$pre.1, alpha = 0.05)) #obtain 95% confidence intervals for the means of a DV
(ci1 <-
    ci.mean(pre.1 ~ gender, data = OBdata)) #obtain 95% confidence intervals for the means of DV for groups on an IV
(ci1 <-
    ci.mean(DV ~ A + B + C, data = data)) #obtain 95% confidence intervals for the means of DV for groups on an IV
plot(ci1) #plot the confidence intervals and means obtained using the ci.means function

#writing your own functions#####
#R is a fully functional programming language
#if a function doesn't exist in the packages available to you, you can write your own
#example: calculating the mean or standard error
#you can use other existing functions in the code :)

######################################################################
#custom function to calculate standard error of the mean
SE <- function(y) {
  print(sd(y) / sqrt(length(y)))
}
#######################################################################

names(data)

SE(OBdata$pre.1) #using the custom function

#verify using the formula
sd(OBdata$pre.1) / sqrt(length(OBdata$pre.1)) #calculation using the formula standard deviation / the square root of the sample size

#verify using other functions

library(Rmisc) #contains the summarySE function
summarySE(data = data,
          measurevar = "mpg",
          na.rm = T) #means, std dev, SE and CI for a DV split into groups on an IV. To get the upper and lower bounds of the CI just add or subtract the CI value from the mean

library(psych) #contains the describe and describeBY functions
describe(data$mpg)