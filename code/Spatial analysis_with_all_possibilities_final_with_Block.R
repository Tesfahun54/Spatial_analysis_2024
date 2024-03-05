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
#########################################
# Import phenotype data 
#######################################
set.seed(123)
setwd("/Users/tas286/Documents/Data_geno_pheno_for_selected_trials")

phen = read.csv(file = "/Users/tas286/Documents/GitHub/Database analysis/Pheno_Comb_YldQt_Val_2014.csv") # load data
head(phen)


dim(phen)
## Change the factor variable to factor
set.seed(123)
phen$geno = as.factor(phen$line_name)
# phen$rowf = as.factor(phen$row)
# phen$colf = as.factor(phen$column)
phen$Block = as.factor(phen$Block)
phen$col = phen$column
phen$row = phen$row
phen$loc = as.factor(phen$trial)
phen$Test_weight = phen$test.weight
phen$Grain_Yield = phen$grain.yield
phen$Plant_height = phen$plant.height
phen$rep = as.factor(phen$replication)

str(phen)
summary(phen)
colnames(phen)
Traits = c("Test_weight", "Plant_height","Grain_Yield")
dim(phen)
length(levels(phen$loc))
length(levels(phen$geno))
summary(phen)
levels(phen$geno)
Env = levels(phen$loc)
##############################################################
# 1. The cross validation of the  without spatial and without marker
##############################################################
str(phen)
###
# With only block row, and block column 
Acc_wom_wosp = tibble()


#spliting the data into trainign and testing set 
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
    # Creating a folder that contain 5 subset with 100 times with a total of 500
    fold5 = caret::createMultiFolds(y = unique(phen$geno), k = 5, times = 5)
    
    
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
      
      #####################
      # make the dsisgn matrix for the blk_rwo and Blk-col
      # Random factor matrix
      idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
      Z.Blk.test <- model.matrix(~idcol - 1)
      rownames(Z.Blk.test) <- test.data$plot
      
       eval(parse(text = paste("ans <- mmer(",Trait,"~1,
                         random=~ Block,
                         rcov=~vsr(units),
                         data= train.data)")))
       #########################
       # blockrow and blockcol effects 
       befall =  as.matrix(ans$U$Block[[Trait]])
       len_b =as.numeric(levels(test.data$Block))
       blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set
       obs.test = test.data[,c("plot",Trait)]
       head(test.data)
       efftest = blkeff 
       r = cbind(env, Trait, predictability = round(cor(efftest[,1],obs.test[,2], use = "pairwise.complete.obs"),3))
       colnames(r)[3] = "predictability"
      Acc_wom_wosp = rbind(Acc_wom_wosp,r) 
      
    }
  }
}
df1 = as.data.frame(Acc_wom_wosp)
head(df1)
df1$method = "Block"
##################################################################
# Cross validation for with spatial analysis but not marker data
 ###################################################################
Acc_wom_wsp = tibble()
    
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
        # Creating a folder that contain 5 subset with 100 times with a total of 500
        fold5 = caret::createMultiFolds(y = unique(phen$geno), k = 5, times = 10)
        
        
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
          
          #####################
          # make the dsisgn matrix for the blk_rwo and Blk-col
          # Random factor matrix
          idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
          Z.Blk.test <- model.matrix(~idcol - 1)
          rownames(Z.Blk.test) <- test.data$plot
         
          
         eval(parse(text = paste("ans1 <- mmer(",Trait,"~1,
                             random=~ Block+ spl2Da(row,col),
                             rcov=~vsr(units),
                             data= train.data)"))) 
         befall =  as.matrix(ans1$U$Block[[Trait]])
         len_b =as.numeric(levels(test.data$Block))
         blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set
          
          # make a plot to observe the spatial effects found by the spl2D()
          W <- with(test.data,spl2Da(row,col)) # 2D spline incidence matrix
          test.data$spatial <- W$Z$`A:all`%*%ans1$U$`A:all`[[Trait]] # 2D spline BLUPs
          
          obs.test = test.data[,c("plot",Trait)]
          efftest = blkeff + test.data$spatial
          
          r = cbind(env, Trait, predictability = round(cor(efftest[,1],obs.test[,2], use = "pairwise.complete.obs" ),3))
          colnames(r)[3] = "predictability"
          Acc_wom_wsp = rbind(Acc_wom_wsp,r) 
          
          
        }
      }
    }
  df2 = as.data.frame(Acc_wom_wsp)
  df2$method = "Block+Spatial"   
###########################################################
# 3. Genomic selection without spatial analysis 
###########################################################
        
###################################
# Estimate the relationship matrix
###################################

##########################
# import the genotyp data 
##########################
###Read the snp file as dosage
snps = read.csv("~/Documents/Data_geno_pheno_for_selected_trials/Imputed_data/Imputed/YldQT_Val_2014_Imputed_Dosage.csv")
head(snps[,1:10])
dim(snps)

# snps[,1] %in% phen$geno
# snps[,1] = gsub(pattern = ".", replacement = "-", x = snps[,1])
 snps = data.frame(snps)
rownames(snps) <- snps[,1]
snps = snps[,-1]
head(snps[,1:10])
dim(snps)
dim(phen)
rownames(snps) %in% unique(phen$geno)
A1 = A.mat(scale(snps,scale = T, center = T)) # The addative relationship matrix
        dim(A1)
        A1 = (1-0.05)*A1 + (0.05)*diag(length(rownames(A1))) # to void singlualrity 
        all(rownames(A1)%in% phen$geno)
        indm = rownames(A1) %in% phen$geno
        A1 = A1[indm, indm]
        indp = phen$geno %in% rownames(A1)
        phen = phen[indp,]
        all(phen$geno %in% rownames(A1))        


 ############################
# 4. with marker and blocrow and blkcol 
#########################################
        
Acc_wm_wosp_bkrc = tibble()
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
            # Creating a folder that contain 5 subset with 100 times with a total of 500
            fold5 = caret::createMultiFolds(y = unique(phen$geno), k = 5, times = 5)
            
            
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
              
              #####################
              # make the dsisgn matrix for the blk_rwo and Blk-col
              # Random factor matrix
              idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
              Z.Blk.test <- model.matrix(~idcol - 1)
              rownames(Z.Blk.test) <- test.data$plot
              
              
              idgtest <- factor(as.character(test.data[,"geno"]), levels = unique(test.data[,"geno"]))
              Z.geno.test <- model.matrix(~idgtest - 1)
              rownames(Z.geno.test) <- test.data$geno
              
              
              eval(parse(text = paste("ans4 <- mmer(",Trait,"~1,
                                 random=~ Block + vsr(geno,Gu = A1),
                                 rcov=~vsr(units),
                                 data= train.data)")))
              
              genoUef = as.matrix(ans4$U$`u:geno`[[Trait]])
              genoUef= as.data.frame(genoUef)
              genoUef$geno = rownames(genoUef)
              test.data$genoeff = NA
              test.data$blockeff = NA
             
              befall =  as.matrix(ans4$U$Block[[Trait]])
              len_b =as.numeric(levels(test.data$Block))
              blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set
              
              test.data$blockeff =  blkeff
              
              for(g in as.vector(test.data$geno)){
                
                test.data[test.data$geno == g,"genoeff"] = genoUef[genoUef$geno == g, "V1"]
              }
              
              test.data$toteff = test.data[,"blockeff"]  + test.data[,"genoeff"]
              
              r = cbind(env, Trait, predictability = round(cor(test.data[,Trait],test.data[,"toteff"], use = "pairwise.complete.obs"),3))
              colnames(r)[3] <- c("predictability")
              Acc_wm_wosp_bkrc = rbind(Acc_wm_wosp_bkrc,r) 
              
              
            }
          }
        }       
  df4 = as.data.frame(Acc_wm_wosp_bkrc)
  df4$method = "Block+Marker"
         
        
        
        
        
        
        
        
###############################################################################
# 5. Spatial analysis with marker and design matrix 
###############################################################################
 Acc_wm_wsp = tibble()       
        
            
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
                      
                      SL = droplevles(SL[-out_ind,])
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
                # Creating a folder that contain 5 subset with 100 times with a total of 500
                fold5 = caret::createMultiFolds(y = unique(phen$geno), k = 5, times = 5)
                
                
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
                  
                  # test.data[,TraitN] = NA # change the grain yield of the training set to NA value
                  
                  mod_dat = rbind(train.data, test.data) # combine the the data set for analysis
                  
                  #####################
                  # make the dsisgn matrix for the blk_rwo and Blk-col
                  # Random factor matrix
                  dim(test.data)
                  idblock <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
                  Z.Blk.test <- model.matrix(~idblock - 1)
                  
                  rownames(Z.Blk.test) <- test.data$plot
            
                  
                  eval(parse(text = paste("ans3<- mmer(",Trait,"~1,
                                     random=~ Block + vsr(geno, Gu = A1) +
                                       spl2Da(row,col),
                                     rcov=~vsr(units),
                                     data= train.data)")))
                  
            
                  befall =  as.matrix(ans3$U$Block[[Trait]])
                  len_b =as.numeric(levels(test.data$Block))
                  blockeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set
                  
                  # make a plot to observe the spatial effects found by the spl2D()
                  W <- with(test.data,spl2Da(row,col)) # 2D spline incidence matrix
                  test.data$spatial <- W$Z$`A:all`%*%ans3$U$`A:all`[[Trait]] # 2D spline BLUPs
                 
                  genoUef = as.matrix(ans3$U$`u:geno`[[Trait]])
                  genoUef = as.matrix(genoUef[order(rownames(genoUef)),])
                  rownames(genoUef) %in% rownames(A1)
                  id = which(rownames(genoUef) %in% unique(test.data$geno))
                  genoUtest = as.matrix(genoUef[id,])
                  
                  rownames(genoUtest) %in% unique(test.data$geno)
                  dim(genoUtest)
                  # geneff = as.matrix(ans3$U$`u:geno`[[Trait]])
                  # geneff = as.matrix(geneff[order(rownames(geneff)),])
                  # rownames(geneff) %in% rownames(A1)
                  # id = which(rownames(geneff) %in% unique(test.data1$geno))
                  # genefftest = as.matrix(geneff[id,])
                  # rownames(genefftest) %in% unique(test.data1$geno)
                  test.data$genoUtest = NA
                  genoUtest = as.data.frame(genoUtest)
                  genoUtest$geno = rownames(genoUtest)
                  ################################################
                  # putting the genetic effect for the indvidual plot
                  #################################################
                  #place the genoU blup in the table 
                  for(g in as.vector(test.data$geno)){
                    
                    test.data[test.data$geno == g,"genoUtest"] = genoUtest[genoUtest$geno == g, "V1"]
                  }
                  ############################
                  # # putting genotypic effec to the data
                  # for(g in as.vector(test.data$geno)){
                  #   
                  #   test.data[test.data$geno == g,"genoeff"] = genefftest[genefftest$geno == g, "V1"]
                  # }
                  # 
                
                  efftest = blockeff + test.data$spatial + test.data$genoUtest
                  test.data$toteff = efftest
                  
                  
                  r = cbind(env,Trait,
                                predictability = round(cor(test.data[,paste(Trait)],
                                                           test.data[,"toteff"],
                                                           use = "pairwise.complete.obs"),3))
                  colnames(r)[3] = "predictability"
                  Acc_wm_wsp = rbind(Acc_wm_wsp,r) 
                  
                  
                }
                
              }
}

df5 = as.data.frame(Acc_wm_wsp)     
df5$method = "Block+Marker+Spatial"  
head(df5)

dfcomb = rbind(df1,df2,df4,df5)
str(dfcomb)
dfcomb$predictability = as.numeric(dfcomb$predictability)
write.csv(x = dfcomb, file = "/Users/tas286/Documents/GitHub/Genomic selection wheat/spatialanaysis_final/sumpred_YLDQt_val_2014_Block.csv")
##########################
# ploting the graph

## stacked box plot 


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
  ggtitle("Predictability of Trial YLDQt_val_2014")

library(gridExtra)
ggarrange(e2)
grid.arrange(e3)



##############################
# EStimate the mean value of the predictablility

df = read.csv(file = "/Users/tas286/Documents/GitHub/Genomic selection wheat/spatialanaysis_final/Crossvalidation_new_block_all_trials.csv")
head(df)
str(df)
library(dplyr)
library(reshape2)
wid1 = dcast(data = df, formula = env + Trait ~ method, fun.aggregate = mean, value.var = "predictability", na.rm = T)

Meancrs$env = factor(Meancrs$env)
Meancrs$Trait = factor(Meancrs$Trait)
Meancrs$Method = factor(Meancrs$Method)
str(Meancrs)
avc = lm(formula = Predictability ~ env + Trait + Trait:Method + Method, data = Meancrs)
summary(avc)
anova(avc)
emmeans(object = avc, specs = c("env"))
