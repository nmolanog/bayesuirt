rm(list=ls())
library(devtools)
library(roxygen2)
library(here)
#setwd("/media/nicolas/DATOS/Rprojects")
create("bayesuirt")
###insert file with functions in R folder
setwd("./bayesuirt")
document()

setwd("..")
install("bayesuirt")
library(bayesuirt)
