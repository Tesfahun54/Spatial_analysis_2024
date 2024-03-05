###############################################################################
# Estimating variance component, heritability , RMSE, MAD and correlation for the whole data
###############################################################################
#############################################
# Data preparation for the analysis
library(tidyverse)
library(sommer)
library(caret)
library(dplyr)
library(lme4)
library(rrBLUP)
#########################################
# Import phenotype data
#######################################

setwd("/Users/tas286/Documents/GitHub/Genomic_selection_wheat/12trials_Marker-data/Imputed_markerdata/")
files_phen = c("ABBmid_2014_Blacksburg","ABBmid_2014_Warsaw", "ABBmid_2014_Lexington",
               "ABBmid_2014_Woodford", "YldQtl-Val_2014_Lincoln","YldQtl-Val_2014_ClayCenter",
               "YldQtl-Val_2014_Mead" , "YldQtl-Val_2014_Sidney" , "HWWpanel_2012_Mead",
               "CSR-Val_2015_Mead", "HWWpanel_2013_Tipton", "HWWpanel_2012_Tipton")
files_snp = c("ABBmid_2014_Blacksburg","ABBmid_2014_Warsaw", "ABBmid_2014_Lexington",
              "ABBmid_2014_Woodford", "YldQtl-Val_2014_Lincoln","YldQtl-Val_2014_ClayCenter",
              "YldQtl-Val_2014_Mead" , "YldQtl-Val_2014_Sidney" , "HWWpanel_2012_Mead",
              "CSR-Val_2015_Mead", "HWWpanel_2013_Tipton", "HWWpanel_2012_Tipton")



c("PC_SNPs_CSR_Val_2014_Mead.csv.txt")
path_phen = "/Users/tas286/Documents/GitHub/data_geno_pheno/"
path_snps = "~/Documents/GitHub/Genomic_selection_wheat/12trials_Marker-data/Imputed_markerdata/"
methods = c("Block", "Block+Spatial", "Block+Marker", "Block+Marker+Spatial")
method_outputs = tibble()

##############


set.seed(123)
setwd("/Users/tas286/Documents/Data_geno_pheno_for_selected_trials")

phen = read.csv(file = "/Users/tas286/Documents/GitHub/Database analysis/Peno_Wide_plotbase_HWWpanel_2012_Tipton.csv") # load data

head(phen)
dim(phen)
## Change the factor variable to factor
set.seed(123)
phen$geno = as.factor(phen$line_name)
# phen$rowf = as.factor(phen$row)
# phen$colf = as.factor(phen$column)
phen$Blk_row = as.factor(phen$Blk_row)
phen$Blk_column = as.factor(phen$Blk_column)
phen$col = phen$column
phen$row = phen$row
phen$loc = as.factor(phen$Trial.Code)
phen$Test_weight = phen$test.weight
phen$Grain_Yield = phen$grain.yield
phen$Plant_height = phen$plant.height
phen$rep = as.factor(phen$replication)
phen$Block = as.factor(phen$Block)
str(phen)
summary(phen)
colnames(phen)
# "Test_weight","Plant_height",
Traits = c( "Grain_Yield")
dim(phen)
length(levels(phen$loc))
length(levels(phen$geno))
summary(phen)
levels(phen$geno)
Env = levels(phen$loc)
levels(phen$Block)

############################
# 1. without marker and without spatial analysis
##########################

# the tibble to store the predictability output from the model
pr_wo_m_wo_sp = tibble()

# The tibble to store the Root Mean Squared Error (RMSQ)  and Mean Absolute Error (MAD)
Accu_wo_m_wo_sp = tibble()

# heritability estimate
H2_Trait_wo_m_wo_sp = tibble()

#Variance component
VC_Trait_wo_m_wo_sp = tibble()

# AIc value
AIC_trait_wo_m_wo_sp = tibble()

BLUPs_wo_m_wo_sp = tibble()

for(env in Env){
  SL <- subset(x = phen, subset = loc == env) # Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait

  ntt = length((TraitN))
  head(SL)

  for(Trait in TraitN){

    #Choosing the method of outlier testing for replicated and unreplicated trials
    if(length(SL$rep)/length(levels(SL$geno)) <= 1){
      # removing outlier using boxplotstat for unreplicated trials
      out_ind <- which(SL[,paste(Trait)] %in% boxplot.stats(SL[,paste(Trait)])$out)

      if(length(out_ind) == 0){
        SL = SL}else{

          SL = SL[-out_ind,]
        }

    }else{
      #removing outlier for replicated trials
      eval(parse(text = paste("outl1 <- lmer(",Trait,"~(1|geno),
               data= SL)")))

      outlier = which(stats::rstudent(outl1) > 3)
      if(length(outlier) == 0){
        SL = SL}else{

          SL = SL[-outlier,]
        }
    }

      eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~ geno + Block,
               rcov=~vsr(units),
               data= SL)")))
      r = length(levels(as.factor(SL$rep)))
      ss = summary(ans)
      vc = data.frame(ss$varcomp)
      rownames(vc) <- c("geno", "Block","residual")

      REP = 1- (mean((diag(ans$PevU$geno[[Trait]])))/(vc["geno","VarComp"]))
      H2cullis = 1 - (mean((diag(ans$PevU$geno[[Trait]])))/(2*vc["geno","VarComp"]))

      h2 = cbind(location = paste(env),Trait = paste(Trait),
                 h2std = round(vc["geno","VarComp"]/(vc["geno", "VarComp"] +
                                                       (vc["residual", "VarComp"])),3),
                 H2 = vpredict(ans,  h~ V1/(V1+V3) ), H2Cullis = H2cullis, Reptability = REP)

      H2_Trait_wo_m_wo_sp = rbind(H2_Trait_wo_m_wo_sp,h2) # store the heritability in the tibble
      vg = vc["geno","VarComp"]
      ve = vc["residual", "VarComp"]
      vph = vg + (ve/r)
      vcomp = cbind(location = paste(env), Trait = paste(Trait),gen.var = round(vg,3),
                    error.var = round(ve,3),   phen.var = round(vph,3),
                    heritability = round(vc["geno","VarComp"]/(vc["geno", "VarComp"] +       (vc["residual", "VarComp"]/r)),3))

      VC_Trait_wo_m_wo_sp = rbind(VC_Trait_wo_m_wo_sp, vcomp)

      ## AIC values of the model
      aic = cbind(location = paste(env), Trait = paste(Trait), AIC = ans$AIC)
      AIC_trait_wo_m_wo_sp = rbind(AIC_trait_wo_m_wo_sp,aic)

      ###################################################################
      ### predict the BLUPs of the test set
      ###################################################################

      ff = fitted(ans)
      fited = ff$dataWithFitted
      Mean =  cbind(location = paste(env), Trait = paste(Trait), fited[,c("geno",paste(Trait,"fitted", sep = "."),Trait)])
      colnames(Mean)[4] = "fitted"
      colnames(Mean)[5] = "observed"
      BLUPs_wo_m_wo_sp = rbind(BLUPs_wo_m_wo_sp, Mean)
      # Estimate the RMSE (Root mean square error) and MAE (mean absolute error) vlaues between the predicted and observed value
      trt = paste(Trait)
      RM = cbind(location = paste(env),Trait = paste(Trait),RMSE = round(RMSE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T) ,3),
                 MAE = round(MAE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T),3))
      # Estimate the predictability among the estiamted and predicted value

      SpearMancor = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],method = "spearman",
                              use = "pairwise.complete.obs"),3)
      preid = cbind(location = paste(env),Trait = paste(Trait),
                    predictability = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],
                                               use = "pairwise.complete.obs"),3),
                    spearman_cor = SpearMancor) # estimate the correlation of the obseved and the predicted value

      Accu_wo_m_wo_sp = rbind(Accu_wo_m_wo_sp,RM) # store the measure of accuracy in tibble


      pr_wo_m_wo_sp = rbind(pr_wo_m_wo_sp,preid) # store model predictability in tibble

    }
  }


#######################
#2. without marker with spatial analysis
###########################################
# the tibble to store the predictability output from the model
pr_wo_m_w_sp = tibble()

# The tibble to store the Root Mean Squared Error (RMSQ)  and Mean Absolute Error (MAD)
Accu_wo_m_w_sp = tibble()

# heritability estimate
H2_Trait_wo_m_w_sp = tibble()

#Variance component
VC_Trait_wo_m_w_sp = tibble()

# AIc value
AIC_trait_wo_m_w_sp = tibble()

#
BLUPs_wo_m_w_sp = tibble()


for(env in Env){
  SL <- subset(x = phen, subset = loc == env) # Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait

  ntt = length((TraitN))
  head(SL)

  for(Trait in TraitN){

    #Choosing the method of outlier testing for replicated and unreplicated trials
    if(length(SL$rep)/length(levels(SL$geno)) <= 1){
      # removing outlier using boxplotstat for unreplicated trials
      out_ind <- which(SL[,paste(Trait)] %in% boxplot.stats(SL[,paste(Trait)])$out)

      if(length(out_ind) == 0){
        SL = SL}else{

          SL = SL[-out_ind,]
        }

    }else{
      #removing outlier for replicated trials
      eval(parse(text = paste("outl1 <- lmer(",Trait,"~(1|geno),
               data= SL)")))

      outlier = which(stats::rstudent(outl1) > 3)
      if(length(outlier) == 0){
        SL = SL}else{

          SL = SL[-outlier,]
        }
    }

    eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~ geno + Block+
               spl2Da(row,col),
               rcov=~vsr(units),
               data= SL)")))
    r = length(levels(as.factor(SL$rep)))
    ss = summary(ans)
    vc = data.frame(ss$varcomp)
    rownames(vc) <- c("geno", "Block","spatial","residual")

    REP = 1- (mean((diag(ans$PevU$geno[[Trait]])))/(vc["geno","VarComp"]))
    H2cullis = 1 - (mean((diag(ans$PevU$geno[[Trait]])))/(2*vc["geno","VarComp"]))

    h2 = cbind(location = paste(env),Trait = paste(Trait),
               h2std = round(vc["geno","VarComp"]/(vc["geno", "VarComp"] +
                                                     (vc["residual", "VarComp"])),3),
               H2 = vpredict(ans,  h~ V1/(V1+V4) ), H2Cullis = H2cullis, Reptability = REP)

    H2_Trait_wo_m_w_sp = rbind(H2_Trait_wo_m_w_sp,h2) # store the heritability in the tibble
    vg = vc["geno","VarComp"]
    ve = vc["residual", "VarComp"]
    vs = vc["spatial", "VarComp"]
    vph = vg + (ve/r)
    vcomp = cbind(location = paste(env), Trait = paste(Trait),gen.var = round(vg,3),
                  error.var = round(ve,3),   phen.var = round(vph,3), sp.var = round(vs,3),
                  heritability = round(vc["geno","VarComp"]/(vc["geno", "VarComp"] +       (vc["residual", "VarComp"]/r)),3))

    VC_Trait_wo_m_w_sp = rbind(VC_Trait_wo_m_w_sp, vcomp)

    ## AIC values of the model
    aic = cbind(location = paste(env), Trait = paste(Trait), AIC = ans$AIC)
    AIC_trait_wo_m_w_sp = rbind(AIC_trait_wo_m_w_sp,aic)

    ###################################################################
    ### predict the BLUPs of the test set
    ###################################################################

    ff = fitted(ans)
    fited = ff$dataWithFitted
    Mean =  cbind(location = paste(env), Trait = paste(Trait), fited[,c("geno",paste(Trait,"fitted", sep = "."),Trait)])
    colnames(Mean)[4] = "fitted"
    colnames(Mean)[5] = "observed"
    BLUPs_wo_m_w_sp = rbind(BLUPs_wo_m_w_sp, Mean)
    # Estimate the RMSE (Root mean square error) and MAE (mean absolute error) vlaues between the predicted and observed value
    trt = paste(Trait)
    RM = cbind(location = paste(env),Trait = paste(Trait),RMSE = round(RMSE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T) ,3),
               MAE = round(MAE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T),3))
    # Estimate the predictability among the estiamted and predicted value

    SpearMancor = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],method = "spearman",
                            use = "pairwise.complete.obs"),3)
    preid = cbind(location = paste(env),Trait = paste(Trait),
                  predictability = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],
                                             use = "pairwise.complete.obs"),3),
                  spearman_cor = SpearMancor) # estimate the correlation of the obseved and the predicted value

    Accu_wo_m_w_sp = rbind(Accu_wo_m_w_sp,RM) # store the measure of accuracy in tibble


    pr_wo_m_w_sp = rbind(pr_wo_m_w_sp,preid) # store model predictability in tibble

  }
}


#################################################
# 3. with marker and without spatial analysis
###############################################
# Import marker data information and estimat Amat
#########################
###Read the snp file as dosage
Pheno_Comb_ABBmid_2014_blacksburg_warsaw.csv
snps = read.csv("~/Documents/Data_geno_pheno_for_selected_trials/Imputed_data/Imputed/HWWpanel_2012_Tipton_Imputed_Dosage.csv")
head(snps[,1:10])

dim(snps)
head(phen)
rownames(snps) <- snps[,1]
snps = snps[,-1]
# snps[,1] %in% phen$geno
# snps[,1] = gsub(pattern = "\\.", replacement = "-", x = rownames(snps))
# snps = data.frame(snps)
# rownames(snps) <- snps[,1]
# snps = snps[,-1]
head(snps[,1:10])
head(phen)
dim(snps)
dim(phen)
A1 = A.mat(scale(snps,scale = T, center = T)) # The addative relationship matrix
dim(A1)
A1 = (1-0.05)*A1 + (0.05)*diag(length(rownames(A1))) # to void singlualrity
A1[1:10,1:10]
all(rownames(A1)%in% phen$geno)
indm = which(rownames(A1) %in% phen$geno)
A1 = A1[indm, indm]
indp = which(phen$geno %in% rownames(A1))
phen = droplevels(phen[indp,])
all(unique(phen$geno) %in% rownames(A1))

# the tibble to store the predictability output from the model
pr_w_m_wo_sp = tibble()

# The tibble to store the Root Mean Squared Error (RMSQ)  and Mean Absolute Error (MAD)
Accu_w_m_wo_sp = tibble()

# heritability estimate
H2_Trait_w_m_wo_sp = tibble()

#Variance component
VC_Trait_w_m_wo_sp = tibble()

# AIc value
AIC_trait_w_m_wo_sp = tibble()

#
BLUPs_w_m_wo_sp = tibble()


for(env in Env){
  SL <- subset(x = phen, subset = loc == env) # Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait

  ntt = length((TraitN))
  head(SL)

  for(Trait in TraitN){

    #Choosing the method of outlier testing for replicated and unreplicated trials
    if(length(SL$rep)/length(levels(SL$geno)) <= 1){
      # removing outlier using boxplotstat for unreplicated trials
      out_ind <- which(SL[,paste(Trait)] %in% boxplot.stats(SL[,paste(Trait)])$out)

      if(length(out_ind) == 0){
        SL = SL}else{

          SL = SL[-out_ind,]
        }

    }else{
      #removing outlier for replicated trials
      eval(parse(text = paste("outl1 <- lmer(",Trait,"~(1|geno),
               data= SL)")))

      outlier = which(stats::rstudent(outl1) > 3)
      if(length(outlier) == 0){
        SL = SL}else{

          SL = SL[-outlier,]
        }
    }


    eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~ Block +
                vsr(geno, Gu = A1),
               rcov=~vsr(units),
               data= SL)")))
    r = length(levels(as.factor(SL$rep)))
    ss = summary(ans)
    vc = data.frame(ss$varcomp)
    rownames(vc) <- c( "Block","u:geno","residual")
    REP = 1- (mean((diag(ans$PevU$`u:geno`[[Trait]])))/(vc["u:geno","VarComp"]))
    H2cullis = 1 - (mean((diag(ans$PevU$`u:geno`[[Trait]])))/(2*vc["u:geno","VarComp"]))

    h2 = cbind(location = paste(env),Trait = paste(Trait),
               h2std = round(vc["u:geno","VarComp"]/(vc["u:geno", "VarComp"] +
                                                     (vc["residual", "VarComp"])),3),
               H2 = vpredict(ans,  h~ V2/(V2+V3) ), H2Cullis = H2cullis, Reptability = REP)


    H2_Trait_w_m_wo_sp = rbind(H2_Trait_w_m_wo_sp,h2) # store the heritability in the tibble
    vg = vc["u:geno","VarComp"]
    ve = vc["residual", "VarComp"]
    vph = vg + (ve/r)
    vcomp = cbind(location = paste(env), Trait = paste(Trait),gen.var = round(vg,3),
                  error.var = round(ve,3),   phen.var = round(vph,3),
                  heritability = round(vc["u:geno","VarComp"]/(vc["u:geno", "VarComp"] +       (vc["residual", "VarComp"]/r)),3))

    VC_Trait_w_m_wo_sp = rbind(VC_Trait_w_m_wo_sp, vcomp)

    ## AIC values of the model
    aic = cbind(location = paste(env), Trait = paste(Trait), AIC = ans$AIC)
    AIC_trait_w_m_wo_sp = rbind(AIC_trait_w_m_wo_sp,aic)

    ###################################################################
    ### predict the BLUPs of the test set
    ###################################################################

    ff = fitted(ans)
    fited = ff$dataWithFitted
    Mean =  cbind(location = paste(env), Trait = paste(Trait), fited[,c("geno",paste(Trait,"fitted", sep = "."),Trait)])
    colnames(Mean)[4] = "fitted"
    colnames(Mean)[5] = "observed"
    BLUPs_w_m_wo_sp = rbind(BLUPs_w_m_wo_sp, Mean)
    # Estimate the RMSE (Root mean square error) and MAE (mean absolute error) vlaues between the predicted and observed value
    trt = paste(Trait)
    RM = cbind(location = paste(env),Trait = paste(Trait),RMSE = round(RMSE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T) ,3),
               MAE = round(MAE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T),3))
    # Estimate the predictability among the estiamted and predicted value

    SpearMancor = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],method = "spearman",
                            use = "pairwise.complete.obs"),3)
    preid = cbind(location = paste(env),Trait = paste(Trait),
                  predictability = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],
                                             use = "pairwise.complete.obs"),3),
                  spearman_cor = SpearMancor) # estimate the correlation of the obseved and the predicted value

    Accu_w_m_wo_sp = rbind(Accu_w_m_wo_sp,RM) # store the measure of accuracy in tibble


    pr_w_m_wo_sp = rbind(pr_w_m_wo_sp,preid) # store model predictability in tibble

  }
}




############################
# 4. with marker and spatial analysis
#########################################
# the tibble to store the predictability output from the model
pr_w_m_w_sp = tibble()

# The tibble to store the Root Mean Squared Error (RMSQ)  and Mean Absolute Error (MAD)
Accu_w_m_w_sp = tibble()

# heritability estimate
H2_Trait_w_m_w_sp = tibble()

#Variance component
VC_Trait_w_m_w_sp = tibble()

# AIc value
AIC_trait_w_m_w_sp = tibble()

#Blups mean
BLUPs_w_m_w_sp = tibble()

for(env in Env){
  SL <- droplevels(subset(x = phen, subset = loc == env)) # Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait

  ntt = length((TraitN))
  head(SL)

  for(Trait in TraitN){

    #Choosing the method of outlier testing for replicated and unreplicated trials
    if(length(SL$rep)/length(levels(SL$geno)) <= 1){
      # removing outlier using boxplotstat for unreplicated trials
      out_ind <- which(SL[,paste(Trait)] %in% boxplot.stats(SL[,paste(Trait)])$out)

      if(length(out_ind) == 0){
        SL = droplevels(SL)}else{

          SL = droplevels(SL[-out_ind,])
        }

    }else{
      #removing outlier for replicated trials
      eval(parse(text = paste("outl1 <- lmer(",Trait,"~(1|geno),
               data= SL)")))

      outlier = which(stats::rstudent(outl1) > 3)
      if(length(outlier) == 0){
        SL = droplevels(SL)}else{

          SL = droplevels(SL[-outlier,])
        }
    }

    eval(parse(text = paste("ans <- mmer(",Trait,"~1,
               random=~  Block+
                vsr(geno, Gu = A1) + spl2Da(row,col),getPEV = T,
               rcov=~vsr(units),
               data= SL)")))
    r = length(levels(as.factor(SL$rep)))
    ss = summary(ans)
    vc = data.frame(ss$varcomp)
    rownames(vc) <- c("Block","u:geno","spatial","residual")

    REP = 1- (mean((diag(ans$PevU$`u:geno`[[Trait]])))/(vc["u:geno","VarComp"]))
    H2cullis = 1 - (mean((diag(ans$PevU$`u:geno`[[Trait]])))/(2*vc["u:geno","VarComp"]))

    h2 = cbind(location = paste(env),Method = "w_m_w_sp", Trait = paste(Trait),
               h2std = round(vc["u:geno","VarComp"]/(vc["u:geno", "VarComp"] +
                                                     (vc["residual", "VarComp"])),3),
               H2 = vpredict(ans,  h~ V3/(V3+V4) ), H2cullis = H2cullis, Reptability = REP)


    H2_Trait_w_m_w_sp = rbind(H2_Trait_w_m_w_sp,h2) # store the heritability in the tibble
    vg = vc["u:geno","VarComp"]
    ve = vc["residual", "VarComp"]
    vs = vc["spatial", "VarComp"]
    vph = vg + (ve/r)
    vcomp = cbind(location = paste(env), Method = "w_m_w_sp", Trait = paste(Trait),gen.var = round(vg,3),
                  error.var = round(ve,3),   phen.var = round(vph,3), sp.var = round(vs,3),
                  heritability = round(vc["u:geno","VarComp"]/(vc["u:geno", "VarComp"] +       (vc["residual", "VarComp"]/r)),3))

    VC_Trait_w_m_w_sp = rbind(VC_Trait_w_m_w_sp, vcomp)

    ## AIC values of the model
    aic = cbind(location = paste(env), Method = "w_m_w_sp",Trait = paste(Trait), AIC = ans$AIC)
    AIC_trait_w_m_w_sp = rbind(AIC_trait_w_m_w_sp,aic)

    ###################################################################
    ### predict the BLUPs of the test set
    ###################################################################

    ff = fitted(ans)
    fited = ff$dataWithFitted
    Mean =  cbind(location = paste(env), Method = "w_m_w_sp",Trait = paste(Trait), fited[,c("geno",paste(Trait,"fitted", sep = "."),Trait)])
    colnames(Mean)[5] = "fitted"
    colnames(Mean)[6] = "observed"
    BLUPs_w_m_w_sp = rbind(BLUPs_w_m_w_sp, Mean)
    # Estimate the RMSE (Root mean square error) and MAE (mean absolute error) vlaues between the predicted and observed value
    trt = paste(Trait)
    RM = cbind(location = paste(env),Method = "w_m_w_sp",Trait = paste(Trait),RMSE = round(RMSE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T) ,3),
               MAE = round(MAE(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")], na.rm = T),3))
    # Estimate the predictability among the estiamted and predicted value

    SpearMancor = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],method = "spearman",
                            use = "pairwise.complete.obs"),3)
    preid = cbind(location = paste(env),Method = "w_m_w_sp",Trait = paste(Trait),
                  predictability = round(cor(fited[,Trait], fited[,paste(Trait,"fitted", sep = ".")],
                                             use = "pairwise.complete.obs"),3),
                  spearman_cor = SpearMancor) # estimate the correlation of the obseved and the predicted value

    Accu_w_m_w_sp = rbind(Accu_w_m_w_sp,RM) # store the measure of accuracy in tibble


    pr_w_m_w_sp = rbind(pr_w_m_w_sp,preid) # store model predictability in tibble

  }
}



####################
# saving the data

# the tibble to store the predictability output from the model
out1 = list(pr_wo_m_wo_sp = pr_wo_m_wo_sp,
            Accu_wo_m_wo_sp =Accu_wo_m_wo_sp,
            H2_Trait_wo_m_wo_sp = H2_Trait_wo_m_wo_sp,
            VC_Trait_wo_m_wo_sp = VC_Trait_wo_m_wo_sp,
            AIC_trait_wo_m_wo_sp = AIC_trait_wo_m_wo_sp,
            BLUPs_wo_m_wo_sp  = BLUPs_wo_m_wo_sp,
            pr_wo_m_w_sp = pr_wo_m_w_sp,
            Accu_wo_m_w_sp = Accu_wo_m_w_sp,
            H2_Trait_wo_m_w_sp = H2_Trait_wo_m_w_sp,
            VC_Trait_wo_m_w_sp = VC_Trait_wo_m_w_sp,
            AIC_trait_wo_m_w_sp = AIC_trait_wo_m_w_sp,
            BLUPs_wo_m_w_sp = BLUPs_wo_m_w_sp,
            pr_w_m_wo_sp = pr_w_m_wo_sp,
            Accu_w_m_wo_sp = Accu_w_m_wo_sp,
            H2_Trait_w_m_wo_sp = H2_Trait_w_m_wo_sp,
            VC_Trait_w_m_wo_sp = VC_Trait_w_m_wo_sp,
            AIC_trait_w_m_wo_sp = AIC_trait_w_m_wo_sp,
            BLUPs_w_m_wo_sp  = BLUPs_w_m_wo_sp,
            pr_w_m_w_sp = pr_w_m_w_sp,
            Accu_w_m_w_sp = Accu_w_m_w_sp,
            H2_Trait_w_m_w_sp = H2_Trait_w_m_w_sp,
            VC_Trait_w_m_w_sp = VC_Trait_w_m_w_sp,
            AIC_trait_w_m_w_sp = AIC_trait_w_m_w_sp,
            BLUPs_w_m_w_sp = BLUPs_w_m_w_sp)

library(openxlsx)
write.xlsx(x = out1,file = "/Users/tas286/Documents/GitHub/Genomic_selection_wheat/spatialanaysis_final/summary_res_h2_vc_All_Trials_combined.xlsx", rowNames = T)



####
# Read from the saved output

library(readxl)
library(tibble)
AIC_comb = tibble()
for(i in c("AIC_trait_wo_m_wo_sp","AIC_trait_wo_m_w_sp","AIC_trait_w_m_wo_sp", "AIC_trait_w_m_w_sp")){
dt <- read_excel(path ="/Users/tas286/Documents/GitHub/Genomic_selection_wheat/spatialanaysis_final/summary_res_h2_vc_All_Trials_combined.xlsx",
           sheet = i)
dt = as.data.frame(dt)
AIC_comb <- rbind(AIC_comb, dt)
}
library(reshape2)
str(AIC_comb)
AIC_comb$AIC <- as.numeric(AIC_comb$AIC)
AIC_comb_wide <- dcast(data = AIC_comb, formula = location + Trait~ Method, fun.aggregate = mean, value.var = "AIC", na.rm = T)
head(AIC_comb_wide)
library(ggplot2)
AIC_comb_wide$location = as.factor(AIC_comb_wide$location)
p1 <- ggplot(data = AIC_comb_wide, aes(x = location, y =  w_m_w_sp -  w_m_wo_sp , fill = location)) +
  geom_bar(stat = "identity") +   theme_bw() + theme( legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "Trial names", y = "Difference (w_m_w_sp, w_m_wo_sp)") + facet_grid(~Trait)
p2 <- ggplot(data = AIC_comb_wide, aes(x = location, y =  w_m_w_sp -   wo_m_w_sp , fill = location)) +
  geom_bar(stat = "identity") +   theme_bw() + theme( legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "Trial names", y = "Difference (w_m_w_sp,  wo_m_w_sp)") + facet_grid(~Trait)

p3 <- ggplot(data = AIC_comb_wide, aes(x = location, y =  w_m_w_sp -   wo_m_wo_sp , fill = location)) +
  geom_bar(stat = "identity") +   theme_bw() + theme( legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "Trial names", y = "Difference (w_m_w_sp,  wo_m_wo_sp)") + facet_grid(~Trait)

p4 <- ggplot(data = AIC_comb_wide, aes(x = location, y =  wo_m_w_sp  -   w_m_wo_sp , fill = location)) +
  geom_bar(stat = "identity") +   theme_bw() + theme( legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "Trial names", y = "Difference (wo_m_w_sp, w_m_wo_sp)") + facet_grid(~Trait)
p5 <- ggplot(data = AIC_comb_wide, aes(x = location, y =  wo_m_w_sp  -   wo_m_wo_sp , fill = location)) +
  geom_bar(stat = "identity") +   theme_bw() + theme( legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "Trial names", y = "Difference (wo_m_w_sp,  wo_m_wo_sp)") + facet_grid(~Trait)

p6 <- ggplot(data = AIC_comb_wide, aes(x = location, y =  w_m_wo_sp  -   wo_m_wo_sp , fill = location)) +
  geom_bar(stat = "identity") +   theme_bw() + theme( legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "Trial names", y = "Difference (w_m_wo_sp,  wo_m_wo_sp)") + facet_grid(~Trait)

library(ggpubr)
# theme(axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.title.x = element_blank() )
ggarrange(p1,
          p2 ,
          p3, ncol = 2,nrow = 2)
ggarrange(p1,p5, ncol = 2, nrow = 1, labels = c("A", "B"), hjust = -3)
head(AIC_comb)
ggplot(AIC_comb, mapping = aes(x = Trait, y = AIC, fill = Method)) + geom_boxplot()+theme_bw()

#boxplot for the difference b/n AIC
method <- c("w_m_w_sp",   "w_m_wo_sp", "wo_m_wo_sp", "wo_m_w_sp")

AIDif1 <- as.data.frame(AIC_comb_wide$w_m_wo_sp - AIC_comb_wide$w_m_w_sp )
AIDif1$Trait <- AIC_comb_wide[,"Trait"]
colnames(AIDif1)[1] = "AIC_Diff"
ggplot(data = AIDif1, mapping = aes(x = Trait, y = AIC_Diff, fill = Trait)) +
  geom_boxplot() + labs(y = "AIC difference (Block+Mrker,Block+Marker+Spatial)")+theme_bw()+
  theme(legend.position = "none")

# AIDif2 <- as.data.frame(AIC_comb_wide$wo_m_w_sp - AIC_comb_wide$w_m_w_sp )
# AIDif2$Trait <- AIC_comb_wide[,"Trait"]
# AIDif3 <- as.data.frame(AIC_comb_wide$w_m_w_sp - AIC_comb_wide$wo_m_wo_sp)
# AIDif3$Trait <- AIC_comb_wide[,"Trait"]
# AIDif4 <- as.data.frame(AIC_comb_wide$w_m_wo_sp - AIC_comb_wide$wo_m_w_sp)
# AIDif4$Trait <- AIC_comb_wide[,"Trait"]
# AIDif5 <- as.data.frame(AIC_comb_wide$w_m_wo_sp - AIC_comb_wide$wo_m_wo_sp)
# AIDif5$Trait <- AIC_comb_wide[,"Trait"]
AIDif6 <- as.data.frame(AIC_comb_wide$wo_m_wo_sp - AIC_comb_wide$wo_m_w_sp)
AIDif6$Trait <- AIC_comb_wide[,"Trait"]
colnames(AIDif6)[1] = "AIC_Diff"

ggplot(data = AIDif6, mapping = aes(x = Trait, y = AIC_Diff, fill = Trait)) +
  geom_boxplot() + labs(y = "AIC difference (Block,BlockSpatial)")+theme_bw()+
  theme(legend.position = "none")


#######################
# bargraph for the heritability
##################################
head(hert)
head(H2_Trait_wo_m_wo_sp)
H2_Trait_wo_m_wo_sp$Method = "H2_Trait_wo_m_wo_sp"
colnames(H2_Trait_wo_m_wo_sp)
H2_Trait_wo_m_w_sp$Method = "H2_Trait_wo_m_w_sp"
colnames(H2_Trait_wo_m_w_sp)
H2_Trait_w_m_wo_sp$Method = "H2_Trait_w_m_wo_sp"
colnames(H2_Trait_w_m_wo_sp)
H2_Trait_w_m_w_sp$Method = "H2_Trait_w_m_w_sp"
colnames(H2_Trait_w_m_w_sp)
H2_Trait_w_m_w_sp = H2_Trait_w_m_w_sp[,c(1,3,4,5,6,7,8,2)]
colnames(H2_Trait_w_m_w_sp)[6] = "H2Cullis"
hert = rbind(H2_Trait_wo_m_wo_sp,H2_Trait_wo_m_w_sp, H2_Trait_w_m_wo_sp,H2_Trait_w_m_w_sp)
hert$Trait = as.factor(hert$Trait)
hert$h2std = as.numeric(hert$h2std)
p1 = ggplot(data = hert, aes(x = Trait, y = h2std, fill = Method)) +
  facet_grid(~location)


p1 + geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
ggtitle("Heritability Trial YLDQt_Val_2014")

unique(hert$Method)
head(hert)
hertBlock <- hert[hert$Method == "H2_Trait_wo_m_wo_sp",c("location", "Trait",  "h2std")]
colnames(hertBlock)[3] = "h2std_block"
hertBlockSP <- hert[hert$Method == "H2_Trait_wo_m_w_sp",c("location", "Trait",  "h2std")]
colnames(hertBlockSP)[3] = "h2std_block_Sp"
hertBlockMrk <- hert[hert$Method == "H2_Trait_w_m_wo_sp" ,c("location", "Trait",  "h2std")]
colnames(hertBlockMrk)[3] = "h2std_block_Mark"
hertBlockMrk_sp <- hert[hert$Method == "H2_Trait_w_m_w_sp"  ,c("location", "Trait",  "h2std")]
colnames(hertBlockMrk_sp)[3] = "h2std_block_Mark_sp"

library(ggplot2)
Hert_comb = cbind(hertBlock,hertBlockSP,hertBlockMrk,hertBlockMrk_sp )
colnames(Hert_comb)
Hert_comb1 = Hert_comb[,c(1,2,3,6,9,12)]
p1 <- ggplot(Hert_comb1, mapping = aes(x = h2std_block, y = h2std_block_Sp)) +
  geom_point(color = "blue")+ facet_grid(~Trait)+ geom_abline(intercept = 0, color = "red", size = 1)+
  theme_bw() + labs(x = "Block", y = "Block+Spatial")
Hert_comb1$h2std_block_Mark
# p2 <- ggplot(Hert_comb1, mapping = aes(x = h2std_block, y =  h2std_block_Mark)) +
#   geom_point(color = "blue")+ facet_grid(~Trait)+ geom_abline(intercept = 0, color = "red", size = 1)+
#   theme_bw() + labs(x = "Block", y = "Block+Marker")

p3 <- ggplot(Hert_comb1, mapping = aes(x = h2std_block_Mark, y = h2std_block_Mark_sp )) +
  geom_point(color = "blue")+ facet_grid(~Trait)+ geom_abline(intercept = 0, color = "red", size = 1)+
  theme_bw() + labs(x = "Block+Marker", y = "Block+Marker+Spatial")


library(ggpubr)
ggarrange(p3,p1 , ncol = 1, nrow = 2, labels = c("A", "B"), hjust = -2 )


#####
# AIC
AIC_trait_wo_m_wo_sp$Method = "AIC_Trait_wo_m_wo_sp"
colnames(AIC_trait_wo_m_wo_sp)
AIC_trait_wo_m_w_sp$Method = "AIC_Trait_wo_m_w_sp"
colnames(AIC_trait_wo_m_w_sp)
AIC_trait_w_m_wo_sp$Method = "AIC_Trait_w_m_wo_sp"
colnames(AIC_trait_w_m_wo_sp)
AIC_trait_w_m_w_sp$Method = "AIC_Trait_w_m_w_sp"
colnames(AIC_trait_w_m_w_sp)

colnames(H2_Trait_w_m_w_sp)[6] = "H2Cullis"
AIC = rbind(AIC_trait_wo_m_wo_sp,AIC_trait_wo_m_w_sp, AIC_trait_w_m_wo_sp,AIC_trait_w_m_w_sp)

library(reshape2)
AIC$AIC = as.numeric(AIC$AIC)
AIC_wide <- dcast(data = AIC, formula = location + Trait ~ Method, fun.aggregate = mean, value.var = "AIC" )

p4 <- ggplot(AIC_wide, mapping = aes(x =AIC_Trait_wo_m_wo_sp, y = AIC_Trait_wo_m_w_sp)) +
  geom_point(color = "blue")+ facet_grid(~Trait)+ geom_abline(intercept = 0, color = "red", size = 1)+
  theme_bw() + labs(x = "Block", y = "Block+Spatial")


p5 <- ggplot(AIC_wide, mapping = aes(x =AIC_Trait_w_m_wo_sp, y = AIC_Trait_w_m_w_sp)) +
  geom_point(color = "blue")+ facet_grid(~Trait)+ geom_abline(intercept = 0, color = "red", size = 1)+
  theme_bw() + labs(x = "Block+Marker", y = "Block+Marker+Spatial")

ggarrange(p4,p5 , ncol = 1, nrow = 2, labels = c("A", "B"), hjust = -2 )




#######################
# ploting the correlation between predicted by differnt methods

head(BLUPs_wo_m_w_sp)

m1 = BLUPs_w_m_w_sp[,c("geno","fitted")]
m2 = BLUPs_wo_m_w_sp[, c("geno","fitted")]
id1 = which(m1$geno %in% m2$geno)
id2 = which(m2$geno %in% m1$geno)
m1 = m1[id1,]
m2 = m2[id2,]

cbind(m1,m2)
cor(m1$fitted,m2$fitted,method = "spearman")
BLUPs_mean$location = as.factor(BLUPs_mean$location)
levels(BLUPs_mean$location)
df1 = BLUPs_mean[BLUPs_mean$location == "ABBmid_2014_Blacksburg",]
df2 = BLUPs_mean[BLUPs_mean$location == "ABBmid_2014_Warsaw",]
library(reshape2)
M1 = dcast(data = df1, formula = geno ~ Method, value.var = "fitted", fun.aggregate = mean )
plot(M1$BLUPs_w_m_w_sp, M1$BLUPs_w_m_wo_sp, col = "blue", pch = 16)
cor(M1$BLUPs_w_m_w_sp, M1$BLUPs_wo_m_w_sp, use = "pairwise.complete.obs")
plot(M1$BLUPs_w_m_w_sp, M1$BLUPs_w_m_wo_sp, col = "blue", pch = 16)

######
# importing data from T3Wheat database
library(devtools)
install_github("TriticeaeToolbox/BrAPI.R")
library(BrAPI)

# Use a known BrAPI Server
conn <- getBrAPIConnection("T3/Wheat")

# Manually set the BrAPI Server host
wheat <- createBrAPIConnection("wheat.triticeaetoolbox.org")


# You only need to do this if you don't want to use version 2 of the BrAPI specification
wheat <- createBrAPIConnection("wheatcap.triticeaetoolbox.org", version="v1")
r = conn$get("germplasm", query=list(germplasmName="JERRY"))
jerry = r$data[[1]]

r = conn$get("observations", query=list(studyDbId=9411))
r$data


####
dt <- read.csv("~/Documents/GitHub/Genomic_selection_wheat/50Trial_data/field_trials_with_layouts_accessions.csv")
head(dt)
unique(dt$trial_name)
Trial_30geno = tibble()
for(i in unique(dt$trial_name)){
  dt1 <- dt[dt$trial_name == i, ]
  geleng <- length(unique(dt1$accession_name))
  if(geleng > 30){
    gl30 <- cbind(Trial = i, No_genotype = geleng)
    Trial_30geno = rbind(Trial_30geno,gl30)
  }
}

dt50 <- read.csv("~/Documents/GitHub/Genomic_selection_wheat/50Trial_data/filtered_complete_dat.csv")
head(dt50)
unique(dt50$studyName)
unique(AIC_comb$location)

