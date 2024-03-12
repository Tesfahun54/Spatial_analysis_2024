###############################################################################
# Function for Estimating variance component, heritability , RMSE, MAD and correlation for the whole data
###############################################################################
#############################################
# Data preparation for the analysis
library(tidyverse)
library(sommer)
library(caret)
library(dplyr)
library(lme4)
library(rrBLUP)
library(tibble)
#########################################
# Import phenotype data
#######################################

source("RCBD_No_spatial.R")
source("RCBD_with_spatial.R")
source("AlphaLattice_without_spatial.R")
source("AlphaLattice_with_spatial.R")
##############

############################
# 1. without marker and without spatial analysis
##########################

spatial_analysis <- function(trait, design1, data, spatial){
  loc <- unique(data$studyName)
  var_comp <- tibble()
  Aic_val <- tibble()
  BLUPs_mean <- tibble()
  for(l in loc){
  #Single location anlaysis '


  loc_data <- droplevels(data[data$studyName == l, ])
  loc_data$blockNumber = as.factor( loc_data$blockNumber )
  loc_data$Block = as.factor( loc_data$Block)
  loc_data$germplasmName = as.factor( loc_data$germplasmName)

  Traits <- trait
  TraitN = colnames(loc_data[Traits])[colSums(is.na(loc_data[Traits])) < 25] # selecting the trait
for(Trait in TraitN ){
 if(design1 == "rcbd" & spatial == "no"){
  output <- RCBD_No_spatial(Trait,loc_data)

 }else if(design1 == "rcbd" & spatial == "yes"){
   output <- RCBD_with_spatial(Trait, loc_data)

  } else if(design1 == "alpha" & spatial == "no"){

    output <- AlphaLattice_without_spatial(Trait,loc_data)

  }else if(design1 == "alpha" & spatial == "yes"){
    output <- AlphaLattice_with_spatial(Trait,loc_data)
}

}

  }

  return(output)
}
