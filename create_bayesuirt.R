rm(list=ls())
library(devtools)
library(roxygen2)
library(here)
#setwd("/media/nicolas/DATOS/Rprojects/bayesuirt")
#create("bayesuirt")
###insert file with functions in R folder
document()

setwd("..")
install("bayesuirt")
library(bayesuirt)
