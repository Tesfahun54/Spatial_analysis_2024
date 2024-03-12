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


##############


set.seed(123)
setwd("/Users/tas286/Documents/Data_geno_pheno_for_selected_trials")

phen = read.csv(file = "data/filtered_complete_dat_edited_as_csv.csv") # load data

summary(phen)

################################
## Change the factor variable to factor

set.seed(123)
factor_names <- c("replicate", "Block", "germplasmName", "studyName")
for(factor in factor_names){
  phen[,factor] = as.factor(phen[,factor])
}


Traits = c( "Grain_Yield", "Plant_height","Test_weight")



###########################
# Disease data transformation using Box-Cox method

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


    try(eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~ germplasmName + Block,
               rcov=~vsr(units),
               data= loc_data)"))), silent = TRUE)

    ################
    # Estimate the variance component, heritability and AIC
    r = length(levels(as.factor(loc_data$rep)))
    ss = summary(ans)
    vc = data.frame(ss$varcomp)
    rownames(vc) <- c("geno", "Block","residual")

    REP = 1- (mean((diag(ans$PevU$geno[[Trait]])))/(vc["geno","VarComp"]))
    H2cullis = 1 - (mean((diag(ans$PevU$geno[[Trait]])))/(2*vc["geno","VarComp"]))

    h2 = cbind(location = paste(env),Trait = paste(Trait),
               h2std = round(vc["geno","VarComp"]/(vc["geno", "VarComp"] +
                                                     (vc["residual", "VarComp"])),3),
               H2 = vpredict(ans,  h~ V1/(V1+V3) ), H2Cullis = H2cullis, Reptability = REP)

    vg = vc["geno","VarComp"]
    ve = vc["residual", "VarComp"]
    vph = vg + (ve/r)
    H2 = round( vg/vph,3)
    vcomp = cbind(location = paste(l), Trait = paste(Trait),gen.var = round(vg,3),
                  error.var = round(ve,3),   phen.var = round(vph,3), heritability = H2)

    var_comp = rbind(var_comp, vcomp)

    ## AIC values of the model
    aic = cbind(location = paste(env), Trait = paste(Trait), AIC = ans$AIC, Design = design1)
    Aic_val = rbind(Aic_val,aic)

    ####
    # BLUPs
    BLUPs <- as.data.frame(ans$U$germplasmName) + ans$Beta$Estimate
    BLUPs$Accession <- rownames(BLUPs)
    colnames(BLUPs) <- c("BLUP", "Accession")
    rownames(BLUPs) <- NULL
    BLUPs <- BLUPs[,c(2,1)]
    BLUPs$Accession <- gsub(pattern = "germplasmName", replacement = "", x = BLUPs$Accession )
    BLUPs <- cbind(location = paste(l), Trait = Trait, BLUPs)

    BLUPs_mean <- rbind(BLUPs_mean, BLUPs)
    output <- list(Aic_val, var_comp, BLUPs_mean)

 }else if(design1 == "rcbd" & spatial == "yes"){
   try(eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~ germplasmName + Block + spl2Da(rowNumber,colNumber) ,
               rcov=~vsr(units),
               data= loc_data)"))), silent = TRUE)

   ################
   # Estimate the variance component, heritability and AIC
   r = length(levels(as.factor(loc_data$rep)))
   ss = summary(ans)
   vc = data.frame(ss$varcomp)
   rownames(vc) <- c("geno", "Block","All", "residual")

   REP = 1- (mean((diag(ans$PevU$geno[[Trait]])))/(vc["geno","VarComp"]))
   H2cullis = 1 - (mean((diag(ans$PevU$geno[[Trait]])))/(2*vc["geno","VarComp"]))

   h2 = cbind(location = paste(env),Trait = paste(Trait),
              h2std = round(vc["geno","VarComp"]/(vc["geno", "VarComp"] +
                                                    (vc["residual", "VarComp"])),3),
              H2 = vpredict(ans,  h~ V1/(V1+V3) ), H2Cullis = H2cullis, Reptability = REP)

   vg = vc["geno","VarComp"]
   ve = vc["residual", "VarComp"]
   vph = vg + (ve/r)
   H2 = round( vg/vph,3)
   vcomp = cbind(location = paste(l), Trait = paste(Trait),gen.var = round(vg,3),
                 error.var = round(ve,3),   phen.var = round(vph,3), heritability = H2)

   var_comp = rbind(var_comp, vcomp)

   ## AIC values of the model
   aic = cbind(location = paste(env), Trait = paste(Trait), AIC = ans$AIC, Design = design1)
   Aic_val = rbind(Aic_val,aic)
   ####
   # BLUPs
   BLUPs <- as.data.frame(ans$U$germplasmName) + ans$Beta$Estimate
   BLUPs$Accession <- rownames(BLUPs)
   colnames(BLUPs) <- c("BLUP", "Accession")
   rownames(BLUPs) <- NULL
   BLUPs <- BLUPs[,c(2,1)]
   BLUPs$Accession <- gsub(pattern = "germplasmName", replacement = "", x = BLUPs$Accession )
   BLUPs <- cbind(location = paste(l), Trait = Trait, BLUPs)

   BLUPs_mean <- rbind(BLUPs_mean, BLUPs)
   output <- list(Aic_val, var_comp, BLUPs_mean)

  } else if(design1 == "alpha" & spatial == "no"){
    loc_data$rep = as.factor( loc_data$replicate)
    loc_data$R = as.factor( loc_data$rowNumber)
    loc_data$C = as.factor( loc_data$colNumber)
    try(eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~ germplasmName + rep + R:rep + C:rep ,
               rcov=~vsr(units),
               data= loc_data)"))), silent = TRUE)

    ################
    # Estimate the variance component, heritability and AIC
    r = length(levels(as.factor(SL$rep)))
    ss = summary(ans)
    vc = data.frame(ss$varcomp)
    rownames(vc) <- c("geno", "rep","R:rep","C:rep","residual")


    vg = vc["geno","VarComp"]
    ve = vc["residual", "VarComp"]
    vph = vg + (ve/r)
    H2 = round( vg/vph,3)
    vcomp = cbind(location = paste(l), Trait = paste(Trait),gen.var = round(vg,3),
                  error.var = round(ve,3),   phen.var = round(vph,3), heritability = H2)

    var_comp = rbind(var_comp, vcomp)

    ## AIC values of the model
    aic = cbind(location = paste(l), Trait = paste(Trait), AIC = ans$AIC, Design = design1)
    Aic_val = rbind(Aic_val,aic)
    ####
    # BLUPs
    BLUPs <- as.data.frame(ans$U$germplasmName) + ans$Beta$Estimate
    BLUPs$Accession <- rownames(BLUPs)
    colnames(BLUPs) <- c("BLUP", "Accession")
    rownames(BLUPs) <- NULL
    BLUPs <- BLUPs[,c(2,1)]
    BLUPs$Accession <- gsub(pattern = "germplasmName", replacement = "", x = BLUPs$Accession )
    BLUPs <- cbind(location = paste(l), Trait = Trait, BLUPs)

    BLUPs_mean <- rbind(BLUPs_mean, BLUPs)
    output <- list(Aic_val, var_comp, BLUPs_mean)

  }else if(design1 == "alpha" & spatial == "yes"){
    loc_data$rep = as.factor( loc_data$replicate)
    loc_data$R = as.factor( loc_data$rowNumber)
    loc_data$C = as.factor( loc_data$colNumber)
    try(eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~ germplasmName + rep + R:rep + C:rep +   spl2Da(rowNumber,colNumber) ,
               rcov=~vsr(units),
               data= loc_data)"))), silent = TRUE)

    ################
    # Estimate the variance component, heritability and AIC
    r = length(levels(as.factor(SL$rep)))
    ss = summary(ans)
    vc = data.frame(ss$varcomp)
    rownames(vc) <- c("geno", "rep","R:rep","C:rep","All","residual")


    vg = vc["geno","VarComp"]
    ve = vc["residual", "VarComp"]
    vph = vg + (ve/r)
    H2 = round( vg/vph,3)
    vcomp = cbind(location = paste(l), Trait = paste(Trait),gen.var = round(vg,3),
                  error.var = round(ve,3),   phen.var = round(vph,3), heritability = H2)

    var_comp = rbind(var_comp, vcomp)

    ## AIC values of the model
    aic = cbind(location = paste(l), Trait = paste(Trait), AIC = ans$AIC, Design = design1)
    Aic_val = rbind(Aic_val,aic)
    ####
    # BLUPs
    BLUPs <- as.data.frame(ans$U$germplasmName) + ans$Beta$Estimate
    BLUPs$Accession <- rownames(BLUPs)
    colnames(BLUPs) <- c("BLUP", "Accession")
    rownames(BLUPs) <- NULL
    BLUPs <- BLUPs[,c(2,1)]
    BLUPs$Accession <- gsub(pattern = "germplasmName", replacement = "", x = BLUPs$Accession )
    BLUPs <- cbind(location = paste(l), Trait = Trait, BLUPs)
    BLUPs_mean <- rbind(BLUPs_mean, BLUPs)
    output <- list(Aic_val, var_comp, BLUPs_mean)
}

}

  }

  return(output)
}

spatial_analysis(trait = c("Grain.test.weight...g.l.CO_321.0001210", "Grain.yield...kg.ha.CO_321.0001218"),spatial = "no",
                 design1 =  "alpha", data = phenex)
unique(phenex$studyName)
phenex <- droplevels(phenex)
summary(phenex)
