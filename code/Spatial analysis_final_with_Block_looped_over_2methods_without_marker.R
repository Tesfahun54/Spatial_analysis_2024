####################################################
# R script for analysis of all possible situations
# 1. without marker and without spatial
# 2. Without marekr but with spatial 
# 3. With marker but without spatial anlaysis
# 4. with marker and with spatial analysis

###############################################
# Data preparation for the analysis 
library(tidyverse)
library(sommer)
library(caret)
library(dplyr)
library(lme4)
library(rrBLUP)
library(readxl)
#########################################
# Import phenotype data 
#######################################
set.seed(123)
setwd("/Users/tas286/Documents/GitHub/Genomic_selection_wheat/50Trial_data/")
files_phen = c("data33_Seleted_includ_Genomic_pred_data.csv")

path_phen = "/Users/tas286/Documents/GitHub/Genomic_selection_wheat/50Trial_data/"
methods = c("Block", "Block+Spatial")
method_outputs = tibble()

#phen = read.csv("/Users/tas286/Documents/GitHub/Genomic_selection_wheat/50Trial_data/data33_Seleted_includ_Genomic_pred_data.csv")
# crosval = read.csv("/Users/tas286/Documents/GitHub/Genomic_selection_wheat/50Trial_data/crossval_block_spat_24_trials.csv")
# unique(crosval$env)
# unique(phen$studyName)
for(file in files_phen){
path_Ph = paste(c(path_phen,file),collapse = "/")
# phen_exl <- read_excel(path_Ph, 
#                       sheet = "all_with_row_col")
# phen = phen_exl
phen = read.csv(file = path_Ph) # load phenotype data

##############################################################
# Edit for the phenotypic and genotypic data for the analysis 
##############################################################

## Change the factor variable to factor
phen$geno = as.factor(phen$germplasmName)
# phen$rowf = as.factor(phen$row)
# phen$colf = as.factor(phen$column)
phen$Block = as.factor(phen$Block)
phen$col = phen$colNumber
phen$row = phen$rowNumber
phen$loc = as.factor(phen$studyName)
phen$Test_weight = phen$Grain.test.weight...g.l.CO_321.0001210
phen$Grain_Yield = phen$Grain.yield...kg.ha.CO_321.0001218
phen$Plant_height = phen$Plant.height...cm.CO_321.0001301
phen$rep = as.factor(phen$replicate)
phen$plot = phen$plotNumber

Trait_name = c("Test_weight", "Plant_height","Grain_Yield")
trait_id = which(Trait_name %in% colnames(phen))
Traits = Trait_name[trait_id]
Env = levels(phen$loc)

########################################
# Fivefold cross validation
#########################################
for(env in Env){
  SL <- droplevels(subset(x = phen, subset = loc == env)) # Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait 
  ntt = length((TraitN))
  head(SL)
  
  if(length(unique(SL$plot)) >= 100){
    SL = SL
  
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
    # Creating a folder that contain 5 subset with 100 times with a total of 500
    fold5 = caret::createMultiFolds(y = unique(SL$geno), k = 5, times = 5)
    
    
    for(i in 1:length(fold5)){
      index = fold5[[i]] # the index of the sample for training set
      #subset the phenotypic data
      train_geno = droplevels(unique(SL$geno)[index])
      train_geno_ind = which(SL$geno %in% train_geno)
      train.data <- droplevels(SL %>%
                                 filter(row_number() %in% train_geno_ind)) # subset the training set
      dim(train.data)
      test.data <- droplevels(SL %>%
                                filter(!row_number() %in% train_geno_ind)) # subset the testing set
      dim(test.data)
      
      #test.data[,TraitN] = NA # change the grain yield of the training set to NA value
      
      mod_dat = rbind(train.data, test.data) # combine the the data set for analysis
      
      ######################## ###########
# Analysis based on the four models 
##########################################
for(method in methods){
  if(method == "Block"){
    
    #####################
    # make the dsisgn matrix for the blk_rwo and Blk-col
    # Random factor matrix
    idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
    Z.Blk.test <- model.matrix(~idcol - 1)
    rownames(Z.Blk.test) <- test.data$plot
    ####
    
    eval(parse(text = paste("ans <- mmer(",Trait,"~1,
                         random=~ Block,
                         rcov=~vsr(units),
                         data= train.data)")))
    
    #predict effects
    #########################
    # blockrow and blockcol effects 
    levels(train.data$Block)
    levels(SL$loc)
    befall =  as.matrix(ans$U$Block[[Trait]])
    len_b =as.numeric(levels(test.data$Block))
    blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set
    obs.test = test.data[,c("plot",Trait)]
    head(test.data)
    efftest = blkeff 
    r = cbind(method, env, Trait, predictability = round(cor(efftest[,1],obs.test[,2], use = "pairwise.complete.obs"),3))
    colnames(r)[4] = "predictability"
    method_outputs = rbind(method_outputs,r) 
    
    }else if(method == "Block+Spatial"){
      
      #Design matrix
      #####################
      # make the dsisgn matrix for the blk_rwo and Blk-col
      # Random factor matrix
      idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
      Z.Blk.test <- model.matrix(~idcol - 1)
      rownames(Z.Blk.test) <- test.data$plot
      
      
      eval(parse(text = paste("ans <- mmer(",Trait,"~1,
                             random=~ Block+ spl2Da(row,col),
                             rcov=~vsr(units),
                             data= train.data)"))) 
      
      #estimate effects
      befall =  as.matrix(ans$U$Block[[Trait]])
      len_b =as.numeric(levels(test.data$Block))
      blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set
      str(test.data)
      # make a plot to observe the spatial effects found by the spl2D()
      W <- with(test.data,spl2Da(row,col)) # 2D spline incidence matrix
      test.data$spatial <- W$Z$`A:all`%*%ans$U$`A:all`[[Trait]] # 2D spline BLUPs
      
      obs.test = test.data[,c("plot",Trait)]
      efftest = blkeff + test.data$spatial
      
      r = cbind(method,env, Trait, predictability = round(cor(efftest[,1],obs.test[,2], use = "pairwise.complete.obs" ),3))
      colnames(r)[4] = "predictability"
      method_outputs = rbind(method_outputs,r) 
    } 
    }
  
      
    }
}
}
  }
}
  str(phen)
  #phen = read.csv("/Users/tas286/Documents/GitHub/Genomic_selection_wheat/50Trial_data/data33_Seleted_includ_Genomic_pred_data.csv")
# checking the number of block in each trial  
  dtt = tibble()
for(i in unique(phen$studyName)){
  sl1 = droplevels(phen[phen$studyName == i, ])
  print(unique(sl1$Block))
  bb = cbind(paste(i),unique(sl1$Block))
  print(bb)
  #dtt = rbind(dtt,bb)
}

  
  
  
unique(method_outputs$env)
write.csv(x = method_outputs, file = "/Users/tas286/Documents/GitHub/Genomic_selection_wheat/50Trial_data/mehtod_outputs_35_trials.csv", row.names = F)
##########################
# ploting the graph

## stacked box plot 

method_outputs$predictability = as.numeric(method_outputs$predictability)
dfcomb = method_outputs
str(method_outputs)
dfcomb$method = factor(x = dfcomb$method, levels = c("Block", "Block+Spatial","Block+Marker", "Block+Marker+Spatial"))
e <- ggplot(dfcomb, aes(x = Trait, y = predictability)) +
  facet_grid(~env) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
e2 <- e + geom_boxplot(
  aes(fill = method),
  position = position_dodge(0.9) 
) + 
  scale_fill_manual(values = c("#999999", "#E69F00", "#100000","#Abc111", "#BCDE2222"))   


e2 + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Predictability of 35 trials")


method_outputs$predictability = as.numeric(method_outputs$predictability)
ggplot(dfcomb, aes(x = Trait, y = predictability, fill = method)) +
  geom_boxplot() + theme_bw()
method_outputs = as.data.frame(method_outputs)
method_outputs$method = as.factor(method_outputs$method)
method_outputs$Trait = as.factor(method_outputs$Trait)
method_outputs$env = as.factor(method_outputs$env)
decom2 = decom1 %>% na.exclude(predictability)
avm = lm(formula = predictability ~ env + method + Trait + method:Trait, data = method_outputs)
anova(avm)
summary(decom1)
emmeans(object = avm, specs = c("env","method"))
library(gridExtra)
ggarrange(e2)
grid.arrange(e3)
Env1 = unique(phen$loc)

for(loc in Env1){
  SL = phen[phen$loc == loc, ]
  print(length(SL$plot))
}
##############################
# EStimate the mean value of the predictablility

df = read.csv(file = "/Users/tas286/Documents/GitHub/Genomic_selection_wheat/50Trial_data/crossval_block_spat_24_trials.csv")
head(df)
df = df[,-1]

crosval_58 = rbind(df,method_outputs)
head(crosval_58)
unique(crosval_58$env)


###############
# GGplot for 58 trial s
ggplot(crosval_58, mapping = aes(x = Trait, y = predictability)) +
  geom_boxplot(aes(fill = method)) + theme_bw()

############
# Analysis of varaince of predictability if there is method, trait or method:trait interaction
crosval_sel = crosval_58 %>% drop_na(predictability)
dim(crosval_sel)
decom1= crosval_sel %>% group_by(env, method, Trait) %>%
  summarise_at(.vars = "predictability", .funs = mean, na.rm = T)

dim(decom1)
avc = lm(formula = predictability ~ env + Trait + Trait:method + method + method:env , data = decom1)
summary(avc)
anova(avc)
library(emmeans)
emmeans(object = avc, specs = "Trait")
emmeans(object = avc, specs = c("Trait","method"))



str(df)
library(dplyr)
library(reshape2)
wid1 = dcast(data = method_outputs, formula = env + Trait + method~ , fun.aggregate = mean, value.var = "predictability", na.rm = T)
d1 = wid1[,1:3]
d1$method = "Block"
d2 = wid1[, c(1,2,4)]
d2$method = "Block+Spatial"
colnames(d1)[3] = "predictability"
colnames(d2)[3] = "predictability"
head(wid1)
dcomb = rbind(d1,d2)
Meancrs$env = factor(Meancrs$env)
Meancrs$Trait = factor(Meancrs$Trait)
Meancrs$Method = factor(Meancrs$Method)
str(Meancrs)
head(crosval)
str(crosval)
decom1= method_outputs %>% group_by(env, method, Trait) %>%
  summarise_at(.vars = "predictability", .funs = median)
dim(decom1)
dcomb1 = dcomb[dcomb$predictability > 0, ]
avc = lm(formula = predictability ~ env + Trait + Trait:method + method + method:env, data = decom1)
summary(avc)
anova(avc)
library(emmeans)
emmeans(object = avc, specs = c("method","Trait"))
mdc = dcomb %>% group_by(env) %>%
  summarise_at(.vars = "predictability", .funs = mean, na.rm = T)
data.frame(mdc)
phen[phen$loc == "Adv_Neo_21",]
