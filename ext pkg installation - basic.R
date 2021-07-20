
# core tools --------------------------------------------------------------
install.packages(
  c("tidyverse", #data manipulation/cleaning/reading/writing: dplyr, tidyr, ggplot2, stringr, purrr, forcats, readr, tibble
    "vroom", #very fast file reading
    "readxl", "xlsx", #read and write excel files
    "devtools", #enables installation (and construction of) packages in development
    "remotes", #install packages from non-CRAN repositories
    "DBI", "RMySQL", "RPostgresSQL", "RSQLite", #import data from relational/SQL databases
    "dbplyr", #translating dplyr syntax to SQL
    "haven", #read and write data from SAS, SPSS, and Stata.
    "readstata13", #read .dta files
    "jsonlite", #read and write JSON tables
    "httr"), #working with HTML connections
  repos = "https://cloud.r-project.org/")

# data transformation & exploratory data analysis -------------------------

#exploratory data analysis
remotes::install_github("bcgov/elucidate", dependencies = TRUE) 

#installing elucidate will also install/update the following dependencies (in
#addition to packages listed above):
# 
# data.table #data manipulation tools optimized for working with large data sets
# DT #interactive tables
# gghalves #half-geoms for ggplot2, used to make raincloud plots with elucidate::plot_raincloud()
# janitor #data cleaning utility functions
# lubridate #dates/times
# patchwork #combining multiple graphs into panels of a complex figure
# plotly #interactive graphics
# trelliscopejs #js facetting for ggplot2 graphs when working with larger datasets

install.packages(c("naniar", #functions for working with/evaluating missing data
                   "survey"), #weighted descriptives for large observational samples
                 repos = "https://cloud.r-project.org/") 

#graphing
install.packages(c(
  "ggsignif","ggExtra", "ggpubr", "gganimate", "gridExtra", "ggcorrplot", "ggsci", "viridis",#ggplot2 extensions
  "extrafont", #import fonts from OS
  "dabestr", #effect estimation plots
  "RColorBrewer", #pre-defined color combinations
  "visualize", #visualize.chisq
  "semPlot", #structural equation model plots
  "gratia", #diagnostic and other plots for generalized additive models: https://www.fromthebottomoftheheap.net/2018/10/23/introducing-gratia/
  "circlize", #chord graphs using chordDiagram()
  "ggmap"), #Download street maps straight from Google maps and use them as a background in your ggplots.
  repos = "https://cloud.r-project.org/")

remotes::install_github("coolbutuseless/ggpattern") #add patterns to ggplots

#tables
install.packages(c(
  "kable", "kableExtra", #make nicer knitr::kables for markdown reports (requires R studio)
  "gt", "formattable", "reactable", #JavaScript/html tables
  "xtable"), 
  repos = "https://cloud.r-project.org/")
#N.B. the xtable function takes an R object (like a data frame) and returns the
#latex or HTML code you need to paste a pretty version of the object into your
#documents.

#Rstudio add-in for building regular expressions: https://www.garrickadenbuie.com/project/regexplain/
remotes::install_github("gadenbuie/regexplain") 

#dplyr meets parallel processing; for datasets with > 10,000,000 rows
#https://github.com/hadley/multidplyr/blob/master/vignettes/multidplyr.md
remotes::install_github("hadley/multidplyr") 


# statistical analysis ----------------------------------------------------

install.packages("psych", repos = "https://cloud.r-project.org/") #factor analysis

#tidyverse style statistical inference and modeling
install.packages("tidymodels", repos = "https://cloud.r-project.org/")

#parametric modelling/inferential stats
install.packages(c(
  "Hmisc", "corrr",#correlations
  "quantreg", #quantile regression
  "car", "afex", "emmeans",#general linear and mixed effects modelling (ANOVA/regression), diagnostics, post-hoc tests
  "heavy", #robust linear (mixed) models for heavy-tailed distributions
  "coxme", #mixed effects Cox models
  "geepack", #generalized estimating equation models
  "lme4", #generalized linear mixed effects models. https://socialsciences.mcmaster.ca/jfox/Courses/soc761/Appendix-Mixed-Models.pdf
  "lmerTest", #F and p values for linear mixed models using Satterthwaite & KR approximations.
  "gamm4", #generalized additive mixed models. Better than mgcv for binary outcome data but lacks residual correlation structures.
  "glmmTMB", #zero-infalted generalized mixed models with capabilities for heteroskedasticity, and residual correlation structures?: https://www.rdocumentation.org/packages/glmmTMB/versions/0.2.3/topics/glmmTMB
  "gvlma", #multiple model diagnostics via gvlma()
  "AICcmodavg", #model selection based on AICc
  "Publish", #confidence intervals for variable means with plotting options ?Publish::ci.mean
  "agricolae", #LSD.test post-hoc comparisons for 3 or fewer comparisons
  "weights"), #weighted chi-square tests, wtd.chi.sq()
  repos = "https://cloud.r-project.org/")

#diagnostics, effect size, and power
install.packages(c(
  "effects", #partial dependence/conditional probabilities and plots for logistic/multinomial regression models
  "performance", #utilities for computing indices of model quality and goodness of fit
  "sjstats", #multiple effect sizes, diagnostics; for linear, generalized and mixed models. https://strengejacke.github.io/sjstats/
  "sjPlot", #diagnostics for linear mixed models
  "epitools", #odds ratios
  "rsq", #full and partial R-squared values for linear models
  "MuMIn", #pseudo-R-squared for mixed models via the r.squaredGLMM() function
  "RLRsim", #restricted likelihood ratio tests for mixed and additive model fit comparisons
  "effects", #partial dependence plots for GLMs models
  "DHARMa", #residual diagnostics for hierarchical (multi-level/mixed) regression models: https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
  "pwr", "Superpower", "simr", #power analysis for r, t, F tests and mixed models
  "rcompanion", "pscl", #McFadden's pseudo-R^2 & likelihood ratio tests for logistic regression models
  "ResourceSelection"), #Hosmer-Lemeshow Goodness of Fit (GOF) Test for logisitic regression
  repos = "https://cloud.r-project.org/")

#non-parametric/"robust" stats
install.packages(c(
  "perm", "permuco", ##permutation tests
  "robust", "statmod", #robust regression (for heteroscedasticity issues)
  "WRS2", #2-way ANOVA with trimmed means or medians; post-hoc tests
  "lmPerm"), #permutation linear models
  repos = "https://cloud.r-project.org/")

#time-series analysis
install.packages(c(
  "zoo","xts", #functions for working with time-series data
  "tsibble", "fable", "feasts", #tidyverse-style time series tools
  "fpp3", "forecast", "forcastHybrid", #time series forecasting, book: https://otexts.com/fpp3/
  "quantmod", "tidyquant"), #Tools for downloading financial data, plotting common charts, and doing technical analysis.
  repos = "https://cloud.r-project.org/")

#Bayesian stats
install.packages(c("BEST", "rstanarm", "brms", "tidybayes"), repos = "https://cloud.r-project.org/")

#r-inla takes quite a while to install so don't bother until you're serious about trying it out
#http://www.r-inla.org
#beginner's guide book: https://www.highstat.com/index.php/beginner-s-guide-to-regression-models-with-spatial-and-temporal-correlation
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# productivity & reporting aids -------------------------------------------
remotes::install_github("ropensci/drake") #workflow reproducibility and organization https://books.ropensci.org/drake/

install.packages("tinytex") #latex
tinytex::install_tinytex()

install.packages(c(
  "beepr", "gmailr", #local (beeping) or remote (e-mail) notifications
  "conflicted", #helps ID and resolve function conflicts between loaded packages
  "doParallel", "foreach", "future.apply", "future",#parallel processing
  "furrr", #parallellized version of purrr
  "sparklyr", "rsparkling", #dplyr backend for spark, supports SQL syntax, and spark/H20 ML functionality for big data https://databricks.com/session/r-and-spark-how-to-analyze-data-using-rstudios-sparklyr-and-h2os-rsparkling-packages
  "shiny", #interactive web apps
  "xaringan", #slides
  "pdftools", #for combining and reading pdf files
  "bookdown", "blogdown"),
  repos = "https://cloud.r-project.org/")

remotes::install_github("gadenbuie/xaringanExtra") #extension for the xaringan package
