#session 1: introduction to R####

#See http://www.statmethods.net/ for many basic functions like graphing and graphical parameter options
#http://www.cookbook-r.com/ the cookbook for R is also a helpful resource
#Shift + left mouse button will activate links in R studio, at least it will in windows...


#For an intro to stats in R, read: Navarro, D. J. (2015) Learning statistics with R: A tutorial for psychology students and other beginners. (Version 0.5) University of Adelaide. Adelaide, Australia
#the companion package, learning statistics with R ("lsr") contains some intuitive and useful basic functions for effect sizes, confidence intervals, Chi-Square tests, etc.

#text after the "#" is used for notes and will be ignored by R. 

#general options/functions####
?function_name # open the help page for the function or other object of interest. 
help(mean)
install.packages("car") #installs specified package, in this case "car" by Sociologist and statistician Dr. John Fox at McMaster.
library(car) #loads the specified package
?car

#change default setup for the project so that packages you use frequently automatically load for that project
getwd()
file.edit("C:/Users/CPH/Desktop/documents/McMaster/research/3xTg-AD Project/experiment 1/R files/.Rprofile")

#Basic operations/calculations#####
5+6
89-28
7000*10
25/5
2^20
sqrt(100)
log(2.3)
log10(100000)

log(1000, base=2)
log10(100)/log10(2)
exp(8)
?exp

#logical operators, mostly useful when selecting subsets of data or programming####
1 == 1 #use double equality symbols to check for equality since "=" is reserved for assignment or value specification
1 == 2

1 != 2 #1 "not equal to 2"
TRUE != FALSE

1 < 2
10 <= 10 #less than or equal to
10 <= 9
10 >= 4 #greater than or equal to
(1 & 3) < 4 #compare both 1 "and(&)" 3 to 4
(1 | 10) < 4 #compare 1 "or(|)" 10 to 4. If at least one comparison is true, R reads the expression as true

#see https://www.statmethods.net/management/operators.html for more examples. 

#data/object types in R####
#class() function will tell you which type or class of data an object is.
?class()
class(4.5) #numeric
class("m") #character
class("4") #if you use quotations, R will read the contents as text
class(TRUE) #or class(T). logical. Mostly used in programming and for some function arguments, e.g. to omit missing values
class(1:3) #integers. The ":" is a shortcut for specifying all numbers between the numbers on each side of it.
class(cbind(1:4, 2:5)) #matrix. Required structure for correlations.
class(iris) #data frame. Most datasets you will work with in R are classified as data frames
class(iris$Species) #factor = categorical variable with multiple levels. The $ symbol is explained below

#convert one data type to another
class(as.matrix(1:3)) #now matrix instead of integer
class(as.character(1:3))
as.character(1:3)

#many datasets are already included in R for educational purposes

?datasets #we will be working primarily with the mt cars and iris datasets for these tutorials
?mtcars
?iris

#Variables####
#pretty much any values/output can be stored/saved in R. This is can enable you to use stored values as arguments for functions.
#The simplest example is creating a variable.

x <- c(1:8) #creates variable 'x' comprised of values between c(). "<-" or "=" is used for assignment
x
x <- 5

rm(x) #delete/remove variable 'x'

x #displays content of variable 'x' 
length(x) #displays the number of elements

class(x) #displays the class of variable 'x', ex. Numeric

x <- cbind(1:8, 9:16) #cbind puts values between commas/parentheses into columns of a matrix. All columns must have the same length
x
class(x)

x[1:4,] #lists values in the 1st 1-4 rows of 'x'. square brackets are used like this to specify [rows, columns] in a data frame, vector, or matrix

x <- rbind(1:8, 9:16) #rbind puts values between commas/parentheses into columns of a matrix. All columns must have the same length
class(x)

x[1:4,] #gives an error now, since we used rbind() and there are only 2 rows
x[1:2,] #lists values in the 1st 1-4 rows of 'x'. square brackets are used like this to specify [rows, columns] in a data frame, vector, or matrix

#transpose a matrix
x <- t(x) #convert rows to columns and vice versa
x[1:4,] #our matrix "x" has been transposed as if we created it using cbind()

#selecting specific values
x[c(2, 5, 6),] #lists the 2nd, 5th, and 6th values of 'x'
x[c(-3, -5, -6),] #lists everthing except the 3rd, 5th, or 6th values of 'x'


names(x) #lists the names of the variables/factors in the data frame. needs to be in data frame format
x <- as.data.frame(x) #converts "x" to a data frame

names(x) <- c("a", "b") #changes the names of the variables in 'df' to 'a', 'b', 'c', and 'd'
names(x)

#Also, the "$" symbol is used in R to specify a variable within a data frame
#this is similar to the "\" in windows explorer for folders or files within a folder

data <- x
length(data$a)

?factor
data$c <- factor(c("a", "b")) #add a grouping variable "c" with levels or groups called "a" and "b". 
#The number of rows must be divisible by the number of factor levels

data$c
levels(data$c) #display the levels in the order that R has them arranged in

#sampling####
?sample #opens the help page for the random sample function
sample(1:40, 6, replace=F) #obtain a random sample 6 numbers without replacement using values ranging between 1 and 40

x <- c(1:30)
x

df <- sample(1:100000, 100, replace=T) #sample with replacement, save it in a vector called "df"

?set.seed

set.seed(seed = 934) #sets criteria for random sampling for variable creation if you want it to be repeatable

random.sample <- rnorm(1000, mean = 100, sd = 1) #create a dataset of random, normally distributed data

#creating a data frame from scratch
y <- c(rnorm(n = 60, mean = 100, sd = 20), rnorm(n = 10, mean = 110, sd = 20)) #creates variable y composed of 60 random scores from a normal distribution with a mean of 100 and sd of 20, along with 10 random scores from a normal distribution with a mean of 110 and sd of 20
g <- factor(rep(seq(1, 7, 1), each = 10), labels = "g", ordered = FALSE) #groups the scores from 'y' into 7 sets (g1,g2,etc) containing 10 scores each
z <- letters[1:5]
data <- as.data.frame(cbind(y, g, z))

class(data)
class(data$y)
data$y <- as.numeric(data$y) #convert to numeric

class(data$z) #check the class of variable "Z"
class(data$g) #check the class of variable "g"

#data structure. The easy way to get structural info on the entire dataset####
str(data)

#renaming a specific variable group/level to make it easier to read
levels(data$g)[levels(data$g)=="1"] <- "level 1" 
levels(data$g)


#changing the order of factor levels (can make graphs look better)
data$z <- factor(data$z, levels = c("a", "c", "b", "d", "e")) #here you specify the order that you want for the "sizes" factor
levels(data$z)

#combine levels of multiple grouping factors into an additional single variable
data$factorGxZ <- interaction(data$g, data$z) #combine factors

levels(data$factorGxZ)

#add a new column to the data frame and name it in the same line####
# install.packages("dpylr") #remove the hashtag/pound symbol at the start of a line to activate code if needed
library(dplyr)

data <- mutate(data, r = y/10) #add a column "C" to the data frame which contains a row-wise proportion of variable A / variable B
data


#get the index number of a specific named column####
names(data)
grep("r", colnames(data)) #useful if you have many columns e.g. 'omics data
 
#exporting data to a .csv file####
write.csv(data, file="practice data.csv") #export to csv file for easy copying
getwd()

#importing data####
#option A - file choice via explorer, easier for single use
imported.data <- read.csv(file.choose(), header=TRUE) #open a browser window to choose the data file, import it, and store it in the generic "data" data frame

#importing data option B - syntax only, faster for repeated use
getwd() # get the working directory
wd <- getwd() #store the working directory in an object called "wd"

list.files(path=wd, pattern="practice data.csv", full.names=TRUE) #list of the files in your working directory

#can copy and paste full path/file_name
data <- read.csv("data.csv", header=TRUE) #import the data file and store it in the generic "data" data frame

#or you can make use of file index number from the list_files object
file_names <- list.files(path=wd, pattern=".csv", full.names=TRUE) #list of the files in your working directory
file_names #print the list of file names, insert full path of desired file into the command below

data <- read.csv(file_names[2], header=TRUE) #import the data file and store it in the generic "data" data frame


#obtaining data subsets####
Adata <- subset(data, z=="a") #extracts the subset of the data frame for group A of variable x and stores it in a new data frame called xA
Bdata <- subset(data, z=="b") #extracts the subset of the data frame for group A of variable x and stores it in a new data frame called xA
Adata
ZAdata <- subset(data, g=="1" & z!="a") #extract data from only the AD genotype control group

g1.g2.data <- subset(data , g=="1" | g=="2") #extract data from only the AD genotype control group

#merging data####
ABdata <- rbind(Adata, Bdata) #combine data frames vertically (by row)

ABdata <- cbind(Adata, Bdata) #combine data frames horizontally (by column)

levels(data$g)
sub1 <- subset(data, g %in% c("level 1", "2")) #extracts multiple subsets of the data from a factor and puts them together in a new data frame "sub 1"

names(data)
data.sub1 <- subset(data[,c(4, 3, 2)]) #extracts only the columns of interest from the original data frame

names(data.sub1)
names(data)
str(data.sub1)

#if you are serious about learning R I highly recommend this introductory course on Data Camp
#https://www.datacamp.com/courses/free-introduction-to-r
