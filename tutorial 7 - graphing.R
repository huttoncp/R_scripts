#session 7: graphing####
#basic graphs and useful graphs from the car package####

#import data####
data <- mtcars #store the mtcars data in a generic dataframe called "data" for convenience

library(car)
OBdata <- OBrienKaiser

#Simple Histogram####
hist(mtcars$mpg)

windows()
bp = barplot(matrix(sample(c(0.5,1),4,replace=T),2,2), col=c("white","grey"), legend=c("A1","A2"), beside=T) #randomly pulls out a number
axis(1, at=colMeans(bp), labels=c("B1", "B2"))

ht.d <- 0.266 #proportion of high testosterone males who are delinquent 
ht.nd <- 1-0.266 #proportion of high testosterone males who are not delinquent
nt.d <- 0.1 #proportion of normal testosterone males who are delinquent 
nt.nd <- 0.9 #proportion of normal testosterone males who are not delinquent
proportions <- rbind(ht.d, ht.nd, nt.d, nt.nd) #combine proportions into a variable (could just enter numbers directly)

windows()
barplot(proportions, main="rates of delinquency in high and normal testosterone males",
        xlab="group", ylab="proportions", col=c("white","grey"), ylim=c(0, 1), 
        names.arg=c("high testosterone\ndelinquent", "high testoseterone\nnon-delinquent",
                    "normal testosterone\ndelinquent", "normal testosterone\nnon-delinquent"),
        beside=TRUE)

#graphing chi-square results####
#install.packages("visualize")

library(visualize)
visualize.chisq(stat = 108.19, df = 1, section = "lower") #insert the observed chi-square and df

#box plots####
names(data)
boxplot(hp ~ gear*cyl, data=data, ylab="horse power", xlab="#cylinders", main="hp by cyl boxplot") #generates a nonparametric box plot using dv and iv from data set "df". 
#Solid line is median, box edges are 25th and 75th percentiles

boxplot(mpg ~ cyl, data=data) 

#scatterplots####
library(car) #contains the scatterplot function and many others - adds marginal boxplots.

names(mtcars)

plot(mpg ~ hp, data = mtcars)
scatterplot(mpg ~ hp, data = mtcars, smooth = F) #scatterplot with regression line. 
#set smooth = T to add lowess smooth line 
#boxplots can be turned off also. See ?scatterplot help page for more options like colours
?scatterplot
windows() #opens a plotting window

names(data)
scatterplot(mpg ~ hp | cyl, data = data, boxplots=T, smooth = F) #split the scatterplot by a grouping factor

library(car)

scatter3d(mpg ~ hp + disp, data = data) #3D scatterplot

names(data)
data.sub <- data[,c(1, 4, 6)] 

pairs(data[,c(1, 4, 6)]) #provides a basic scatterplot matrix

scatterplotMatrix(~ mpg+hp+drat+wt|cyl, data=mtcars, main="Three Cylinder Options") #car scatterplot matrix
?scatterplotMatrix

#3D Scatterplot with Coloring and Vertical Lines and Regression Plane #### 
library(scatterplot3d) 
?scatterplot3d

windows()
s3d <-scatterplot3d(y = data$brn.wt.bdy.wt.ratio, z = data$cells.per.DG, x = data$spleen.wt, pch=16, highlight.3d=TRUE,
                    ylab="brain weight/body weight ratio", zlab = "DCX+ cells per DG", xlab="spleen weight (g)",
                    type="h", angle=45, cex.symbols = 2, cex.lab = 1.25, cex.axis = 1.25) #, angle=60
fit <- lm(data$cells.per.DG ~ data$spleen.wt + data$brn.wt.bdy.wt.ratio)

s3d$plane3d(fit)
?scatterplot3d

#Pie Chart with Percentages####
slices <- c(10, 12, 4, 16, 8) 
lbls <- c("US", "UK", "Australia", "Germany", "France")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

#3D pie chart####
library(plotrix)
slices <- c(10, 12, 4, 16, 8) 
lbls <- c("US", "UK", "Australia", "Germany", "France")
pie3D(slices,labels=lbls,explode=0.1,
      main="Pie Chart of Countries ")

#3D pie chart with percentages####
x <- c(32, 12, 30, 45)
labels <- c("California", "Paris", "Moscow", "Mumbai")
pct <- round(x/sum(x)*100)
lbls <- paste(labels, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
windows()
pie3D(x,labels=lbls,explode=0.1,main="city_pie_chart")
legend("topright", c("California", "Paris", "Moscow", "Mumbai"), cex=0.8,fill=rainbow(length(x)))

###################################################################################

#linear model diagnostic plots####
library(car)
lm.1 <- lm(mpg ~ cyl, data=data) #specify the linear model structure, categorical predictors only

qqPlot(lm.1, distribution="norm", main="QQ Plot", ylab="Studentized Residuals", col="black", col.lines="green", envelope=.95) #graphical estimate of normality with 95% confidence envelope

influenceIndexPlot(lm.1, id.n=6)

windows()
plot(lm.1, which=4) #plot of cook's D to identify potential influencial observations
(cutoff <- 4/(summary(lm.1)$df[2]-1)) #cutoff value for Cook's D is 4/(N-k-1) or 4/(denominator total df - 1)
abline(h=cutoff, col="red", lty=2) #add the cut-off line

#########################################################################
#....ggplot2 graphing............................####
#helpful website 1: http://www.cookbook-r.com/Graphs/
#helpful website 2: contains detailed option info: http://ggplot2.tidyverse.org/reference/
#helpful website 3: https://plot.ly/ggplot2/
#helpful website 4: https://www.r-graph-gallery.com/portfolio/ggplot2-package/
#############################################################################

library(ggplot2)
library(Rmisc) #contains the summarySE and summarySEwithin functions

#density plot####
names(data)
class(data$cyl)
data$cyl <- as.factor(data$cyl)
  
ggplot(data, aes(x=hp, fill=cyl)) + 
  geom_density(alpha=.3)


#barplot####
install.packages("ggplot2")
library(ggplot2)

windows()
data$cyl <- as.factor(data$cyl)
bp<- ggplot(data, aes(x=cyl, y=hp, fill=cyl))+
  geom_bar(width = 1, stat = "identity")
bp


#stacked barplot####
mtcars$cyl <- as.factor(mtcars$cyl)
library(ggplot2)
bp<- ggplot(mtcars, aes(x="", y=hp, fill=cyl))+
  geom_bar(width = 1, stat="identity") + coord_flip()
bp
?geom_bar
#A pie chart = stacked bar chart with polar coordinates####
library(scales)
pie <- ggplot(mtcars, aes(x = factor(1), fill = factor(cyl))) +
  geom_bar(width = 1) +
  coord_polar(theta = "y")
pie

#histogram####
ggplot(data, aes(x=mpg)) +
  geom_histogram(binwidth=.5, colour="black", fill="red")+
  theme_bw()

#Histogram overlaid with kernel density curve####
ggplot(data, aes(x=mpg)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

# Overlaid histograms with means####
# get the mean of each group
data.sum1 <- summarySE(data=data, measurevar="mpg", 
                       groupvars=c("cyl"), na.rm=T) # can also specify removal of NA values if missing data
data.sum1 #view the summary of the data

ggplot(data, aes(x = mpg, fill = cyl)) +
  geom_histogram(binwidth=.5, alpha =.5, position="identity") +
  geom_vline(data = data.sum1, aes(xintercept=mpg,  colour=cyl),
             linetype="dashed", size=1)

# Density plots with means####
ggplot(data, aes(x = mpg, fill = cyl)) +
  geom_density(alpha=.3) +
  geom_vline(data = data.sum1, aes(xintercept=mpg,  colour=cyl),
             linetype="dashed", size=1)

#A basic box plot with the conditions colored####
ggplot(data, aes(x=cyl, y=hp, fill=cyl)) + geom_boxplot() + theme_bw()

#The above adds a redundant legend. With the legend removed
ggplot(data, aes(x=cyl, y=hp, fill=cyl)) + geom_boxplot() +
  guides(fill=FALSE)

#add a diamond at the mean
ggplot(data, aes(x=cyl, y=hp, fill=cyl)) + geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)

#x-y scatterplot with 2 continuous variables and a linear regression line. Useful for visualizing correlations.####
names(data)

library(ggplot2)
fig1 <- ggplot(data) + 
  geom_point(aes(x=hp, y=mpg), shape=8, size=8) + #can customize shape and size of points. See http://www.cookbook-r.com/ for options
  xlab("horsepower") +
  ylab("miles per gallon") + 
  geom_smooth(method=lm, aes(x=hp, y=mpg),se=F) + #if you set SE=T then you will also get a SE confidence envelope around your regression line
  theme_bw()
  windows()
fig1

library(ggExtra) #can add marginal density plots or histograms to a ggplot scatterplot
# Marginal density plot
ggMarginal(fig1 + theme_gray())

# Marginal histogram
ggMarginal(fig1 + theme_gray(), type = "histogram",
           fill = "steelblue", col = "darkblue")

#x-y scatterplot for one DV and different levels of a grouping factor (IV). Useful for visualizing raw data for each group.####
names(data)
library(ggplot2)

fig1 <- ggplot(OBdata) + 
  geom_point(aes(x=treatment, y=post.4), shape=2, size=8) + #can customize shape, size and colour of points. See http://www.cookbook-r.com/ for options
  xlab("treatment group") +
  ylab("post test time 4 score") +
  windows()
fig1

#scatterplot with color and shape by categorical factor and are semi-transparent####
data <- mtcars
str(data)

data$cyl <- as.factor(cyl)
str(data)
levels(data$cyl)

ggplot(data, aes(x=mpg, y=hp, color = cyl, shape = cyl)) +
  geom_point(size=6, alpha=0.6) + #alpha arg sets transparency
  geom_smooth(method=lm, aes(x=mpg, y=hp),se=F) + #if you set SE=T then you will also get a SE confidence envelope around your regression line
  ylab("horsepower") + xlab("miles per gallon") +
  scale_color_discrete(name = "# engine cylinders")+ #color legend title. Add arguments breaks = c("") and levels = c("") to relabel the levels
  scale_shape_discrete(name  ="# engine cylinders") #shape legend title. Add arguments breaks = c("") and levels = c("") to relabel the levels

#correlogram####
devtools::install_github("kassambara/ggcorrplot")
library(ggplot2)
library(ggcorrplot)

#1st get correlation matrix
data(mtcars)
corr <- round(cor(mtcars), 1)

#generate the correlogram
ggcorrplot(corr, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Correlogram of mtcars", 
           ggtheme=theme_bw)

#bargraph of group means +/- SE####

names(data)
library(Rmisc)
data.sum1 <- summarySE(data=data, measurevar="mpg", 
                       groupvars=c("cyl"), na.rm=T) # can also specify removal of NA values if missing data
data.sum1 #view the summary of the data

?geom_bar
bargraph.1 <- ggplot(data.sum1, aes(x=cyl, y=mpg, fill=cyl)) + #insert x and y variable names
  geom_bar(position=position_dodge(),colour="black",stat="identity", size=.3, width=0.75) +
  geom_errorbar(aes(ymin=mpg-se, ymax=mpg+se), #insert y variable name
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  ylab("miles per gallon") +
  xlab("number of cylinders")+
  ggtitle("fuel efficiency of US cars as a function of engine cylinders", subtitle="subtitle goes here")+
  scale_fill_hue()+  #each bar will be a different colour
  ylim(0,50)+ #activate to specify the y-axis limits
  theme_bw()

windows()

bargraph.1

#bargraph with SE error bars using ObrienKaiser data from car package as an example####
names(OBdata) #find out what your DVs and IVs are called (case-sensitive) in the datafile
data.sum <- summarySE(OBdata, measurevar="post.5", groupvars=c("treatment","gender"))
data.sum

library(ggplot2)
levels(OBdata$gender)
ggplot(data.sum, aes(x=treatment, y=post.5, fill=gender)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=post.5-se, ymax=post.5+se), #specify your y-variable for error bars based on SE
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  xlab("treatment") + #x-axis label
  ylab("post-test hour 5 score") + #y-axis label
  scale_fill_hue(name="gender", # Legend label, use darker colors
                 breaks=c("F", "M"),
                 labels=c("female", "male")) + #label the z-factor groups as identified using levels() function
  ggtitle("The Effect of treatment and gender on post-test hour 5 score") + #specify the title you want here
  ylim(0,10)+ #activate to specify the y-axis limits
  theme_bw()


#bargraph with custom-specified colours, use scale_fill_manual instead of scale_fill_hue####
data.sum1 <- summarySE(data=data, measurevar="mpg", 
                       groupvars=c("cyl"), na.rm=T) # can also specify removal of NA values if missing data
data.sum1 #view the summary of the data

levels(data$cyl)

bargraph.1 <- ggplot(data.sum1, aes(x=cyl, y=mpg, fill=cyl)) + #insert x and y variable names
  geom_bar(position=position_dodge(),colour="black",stat="identity", size=.3, width=0.75) +
  geom_errorbar(aes(ymin=mpg-se, ymax=mpg+se), #insert y variable name
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  ylab("miles per gallon") +
  xlab("number of cylinders")+
  ggtitle("fuel efficiency of US cars as a function of engine cylinders")+
  scale_fill_manual(name="# of cylinders", # Legend label, use darker colors	
                    breaks=c("4", "6", "8"),	#specify the groups as ordered and coded by R (use levels function)
                    labels=c("Four", "Six", "Eight"),	#specify what you want the labels to appear as
                    values=c("brown", "blue", "green")) + #for colour options see http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/	
  #ylim(0,100)+ #activate to specify the y-axis limits
  theme_bw()

windows()


bargraph.1

#combine multiple graphs (e.g. scatterplots)####
library(ggplot2)
data4 <- subset(data, cyl=="4", na.rm=TRUE) #extracts the subset of the dataframe for group A of variable x and stores it in a new data frame called xA
data8 <- subset(data, cyl=="8", na.rm=TRUE)


fig1a <- ggplot(data4) + 
  geom_point(aes(x=hp, y=mpg), shape=8, size=8) + #can customize shape and size of points. See http://www.cookbook-r.com/ for options
  xlab("horsepower") +
  ylab("miles per gallon") + 
  ggtitle("4 cylinder cars")+
  geom_smooth(method=lm, aes(x=hp, y=mpg),se=T) + #if you set SE=T then you will also get a SE confidence envelope around your regression line
  windows()
fig1a

fig1b <- ggplot(data8) + 
  geom_point(aes(x=hp, y=mpg), shape=8, size=8) + #can customize shape and size of points. See http://www.cookbook-r.com/ for options
  xlab("horsepower") +
  ylab("miles per gallon") + 
  ggtitle("8 cylinder cars")+
  geom_smooth(method=lm, aes(x=hp, y=mpg),se=T) + #if you set SE=T then you will also get a SE confidence envelope around your regression line
  windows()
fig1b


library(gridExtra)

grid.arrange(fig1a, fig1b, ncol=2)  #compile both plots in the same figure

#linegraph with means +/- SE for repeated measures data####

#1st convert a wide-format data set to long format if needed
library(car)
OBdata <- OBrienKaiser
names(OBdata)

OBdata$ID <- c(1:16) #add a subject ID variable if there isn't one

library(reshape2)
OBdata.long <- melt(OBdata, id.vars=c("ID", "treatment", "gender"), 
                    measure.vars = c("pre.5", "post.5", "fup.5"), variable.name = c("trial"), value.name=c("score"))

#summarize data using the summarySEwithin function 
library(Rmisc)
data.sum <- summarySEwithin(data=OBdata.long, measurevar="score", 
                              betweenvars=c("treatment"), withinvars= c("trial"),
                              idvar= c("ID"), na.rm=T) # can also specify removal of NA values if missing data.sum
data.sum

levels(OBdata$treatment)

fig1 <- ggplot(data.sum, aes(x=trial, y=score, group=treatment, colour=treatment, shape=treatment)) +
  geom_errorbar(aes(ymin=score-se, ymax=score+se), 
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.1))  +
  geom_line(size=1.5, aes(linetype = treatment)) + 
  ylab("measurement score") + geom_point(cex=3) +
  scale_x_discrete(name="trial",
                   breaks=c("pre.5", "post.5", "fup.5"),
                   labels=c("pre-test hr 5", "post-test hr 5", "follow-up hr 5")) +
  scale_shape_discrete(name="Treatment",    # Legend label, use darker colors
                       breaks=c("control", "A", "B"), #what they are to R
                       labels=c("ctrl", "low dose", "high dose"))+
  scale_linetype_discrete(name="Treatment",    # Legend label, use darker colors
                          breaks=c("control", "A", "B"), #what they are to R
                          labels=c("ctrl", "low dose", "high dose"))+
  scale_colour_hue(name="Treatment",    # Legend label, use darker colors
                   breaks=c("control", "A", "B"), #what they are to R
                   labels=c("ctrl", "low dose", "high dose"),#what you want them to be
                   l=40) + # Use darker colors, lightness=40 
   ggtitle("main title", subtitle="subtitle")+
  #ylim(0, 120)+ #activate if you want to specify the y-axis limits
  #theme(legend.position="none")+ #activate to remove the legend
  theme_bw(base_size = 16)
  
windows()
fig1

#linegraph with custom colours####
fig1 <- ggplot(data.sum, aes(x=trial, y=score, group=treatment, colour=treatment, shape=treatment)) +
  geom_errorbar(aes(ymin=score-se, ymax=score+se), 
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.1))  +
  geom_line(size=1.5, aes(linetype = treatment)) + 
  ylab("measurement score") + geom_point(cex=3, fill="white") +
  scale_x_discrete(name="trial",
                   breaks=c("pre.5", "post.5", "fup.5"),
                   labels=c("pre-test hr 5", "post-test hr 5", "follow-up hr 5")) +
  scale_shape_discrete(name="Treatment",    # Legend label, use darker colors
                       breaks=c("control", "A", "B"), #what they are to R
                       labels=c("ctrl", "low dose", "high dose"))+
  scale_linetype_discrete(name="Treatment",    # Legend label, use darker colors
                          breaks=c("control", "A", "B"), #what they are to R
                          labels=c("ctrl", "low dose", "high dose"))+
  scale_color_manual(name="Treatment", breaks=c("control", "A", "B"), #what they are to R
                     labels=c("ctrl", "low dose", "high dose"),
                     values=c("red", "blue", "purple"))+
  
  ggtitle("main title", subtitle="subtitle")+
  #ylim(0, 120)+ #activate if you want to specify the y-axis limits
  #theme(legend.position="none")+ #activate to remove the legend
  theme_bw(base_size = 16)

windows()

fig1

#INTERACTIVE GRAPHICS####

library(dplyr)
library(ggplot2)
library(plotly)

glimpse(mtcars)

plot1 <- mtcars %>% ggplot(aes(x = hp, y = mpg, colour = factor(cyl))) +
                             geom_point(position = "jitter")
plot1

ggplotly(plot1)


plot1 <- mtcars %>% ggplot(aes(x = mpg, fill = factor(cyl))) +
  geom_density(alpha = 0.6)
plot1

ggplotly(plot1)

#significance indicators####
library(ggsignif)
windows()
restraint_data %>% 
  na.omit() %>% 
  group_by(group) %>% 
  summarize(mean = mean(restraint_score), 
            se = sd(restraint_score, na.rm = T)/sqrt(length(restraint_score))) %>% 
  ggplot(aes(x = group, y = mean, color = group)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), 
                position = "dodge", width = 0.25, size = 2) +
  geom_point(cex = 7) +
  labs(y = "Mean restraint score by group", x = "Group")+
  ylim(0, 2) +
  scale_x_discrete(name = "Group", breaks = c("SH", "TBI"), labels = c("Sham", "TBI")) +
  scale_color_manual(name = "Group", 
                     breaks = c("SH", "TBI"),
                     values = c("blue", "red"),
                     labels = c("Sham", "TBI")) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(face = "bold")) +
  geom_signif(comparisons = list(c("SH", "TBI")),  annotations="* p < 0.05", 
              color='black', size = 2, textsize = 5,
              y_position = 1.75, map_signif_level=TRUE, na.rm = T)

#purrr::map plotting####
#boxplots
EPM_data %>% 
  select_if(is.numeric) %>%
  map(~boxplot(. ~ EPM_data$group, data = EPM_data))


#purrr::map ggplotting####
var_labs <- c("open arm entries", "closed arm entries", "total arm entries", 
              "risk assessments", "open arm time (s)", "# of fecal boli")

plots <- EPM_data %>% select(open.arm.entries:fecal.boli) %>% 
  select_if(is.numeric) %>%
  map2(var_labs,
       ~ggplot(EPM_data, mapping = aes(y = .x, x = group)) +
         geom_boxplot() +
         geom_jitter(cex = 3, aes(y = .x, x = group, color = convulsions, shape = sex)) + 
         ylab(.y) +
         theme_bw(base_size = 14)
  )

#combine plots

library(Rmisc)
windows()
multiplot(plotlist = plots, cols=3)


