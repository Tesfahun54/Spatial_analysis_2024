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

files_phen = read.csv("data/phenotype_35_Trials_second_times.csv")
unique(files_phen$)

files_snp = c("SNPs_GP_502_2022.tsv",
              "SNPs_SW-AMPanel_2014_90k.tsv", "SNPs_SW-AMPanel_2014_9k.tsv", "SNPs_TP_2018.tsv" ,
              "SNPs_YLDqt_2014.tsv","SNPs_YLDQtl-Val_2015.tsv")

library(stringr)
unique(files_phen$studyName)
id <- which(str_detect(string = files_phen$studyName, pattern = "SW-AMPanel_2014"))
Phen_sel <- files_phen[id,]

library(readr)
Snps <- read_tsv(file = "data/SNPs_SW-AMPanel_2014_9k.tsv")
head(Snps[,1:10])
Snps = as.data.frame(Snps)
rownames(Snps) <- Snps$Marker
Snps = Snps[,-1]
Snpst <- t(Snps)
head(Snpst[,1:10])
dim(Snpst)
Snp_NA = c()
for(i in 1:ncol(Snpst)){
  if(length(which(is.na(Snpst[,i]))) > 0.2*nrow(Snpst)){
    Snp_NA <- c(Snp_NA,i)
  }
}

if(length(Snp_NA) == 0){
  Snpst_wo_NA <- Snpst
}else{
  Snpst_wo_NA <- Snpst[,-Snp_NA]
}


head(Snpst_wo_NA[,1:10])
library(genomicMateSelectR)
Snpst_MAF <- maf_filter(M = Snpst_wo_NA, thresh = 0.05)
dim(Snpst_MAF)


A1a = A.mat(scale(Snpst_MAF,scale = T, center = T)) # The addative relationship matrix
dim(A1a)
A1b = A.mat(scale(Snpst_MAF,scale = T, center = T)) # The addative relationship matrix
dim(A1b)
K2 <- list(A1a,A1b)
A1comb  <- CovComb(Klist = K2)
# library(CovCombR)
A1 = A1comb


A1 = (1-0.05)*A1 + (0.05)*diag(length(rownames(A1))) # to void singlualrity
A1[1:10, 1:10]

idk <- which(Phen_sel$germplasmName %in% rownames(A1))

Phen_sel1 <- Phen_sel[idk,]
dim(Phen_sel1)
all( unique(Phen_sel1$germplasmName) %in% unique(rownames(A1)))
all(unique(rownames(A1) %in% unique(Phen_sel1$germplasmName)))
##

methods = c("Block", "Block+Spatial", "Block+Marker", "Block+Marker+Spatial")

Traits <- c("Grain.yield...kg.ha.CO_321.0001218", "Plant.height...cm.CO_321.0001301",
            "Grain.test.weight...g.l.CO_321.0001210")

library(MASS)
x <- Phen_sel1$Stripe.rust.severity.....CO_321.0001396[which(!is.na(Phen_sel1$Stripe.rust.severity.....CO_321.0001396))]
lm1 <- lm( Stripe.rust.severity.....CO_321.0001396~geno,data = Phen_sel1 )

boxcox(object = lm1$coefficients)

Env <- unique(Phen_sel1$studyName)
summary(!Phen_sel1$Stripe.rust.severity.....CO_321.0001396)

Phen_sel1$loc = as.factor(Phen_sel1$studyName)
Phen_sel1$geno <- as.factor(Phen_sel1$germplasmName)
Phen_sel1$Block = as.factor(Phen_sel1$Block)
Phen_sel1$row <- Phen_sel1$rowNumber
Phen_sel1$col <- Phen_sel1$colNumber
Phen_sel1$plot <- Phen_sel1$plotNumber
Phen <- Phen_sel1

method_outputs = tibble()
    ########################################
    # Fivefold cross validation
    #########################################
    for(env in Env){
      SL <- subset(x = Phen, subset = studyName == env) # Subseting the data for each env
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
        fold5 = caret::createMultiFolds(y = unique(Phen$geno), k = 5, times = 5)


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

              # make a plot to observe the spatial effects found by the spl2D()
              W <- with(test.data,spl2Da(row,col)) # 2D spline incidence matrix
              test.data$spatial <- W$Z$`A:all`%*%ans$U$`A:all`[[Trait]] # 2D spline BLUPs

              obs.test = test.data[,c("plot",Trait)]
              efftest = blkeff + test.data$spatial

              r = cbind(method,env, Trait, predictability = round(cor(efftest[,1],obs.test[,2], use = "pairwise.complete.obs" ),3))
              colnames(r)[4] = "predictability"
              method_outputs = rbind(method_outputs,r)

            }else if(method == "Block+Marker"){
              #Design effect
              #####################
              # make the dsisgn matrix for the blk_rwo and Blk-col
              # Random factor matrix
              idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
              Z.Blk.test <- model.matrix(~idcol - 1)
              rownames(Z.Blk.test) <- test.data$plot


              idgtest <- factor(as.character(test.data[,"geno"]), levels = unique(test.data[,"geno"]))
              Z.geno.test <- model.matrix(~idgtest - 1)
              rownames(Z.geno.test) <- test.data$geno


              eval(parse(text = paste("ans <- mmer(",Trait,"~1,
                                 random=~ Block + vsr(geno,Gu = A1),
                                 rcov=~vsr(units),tolParInv = 0.0019,
                                 data= train.data)")))

              #estimationof effects
              genoUef = as.matrix(ans$U$`u:geno`[[Trait]])
              genoUef= as.data.frame(genoUef)
              genoUef$geno = rownames(genoUef)
              test.data$genoeff = NA
              test.data$blockeff = NA

              befall =  as.matrix(ans$U$Block[[Trait]])
              len_b =as.numeric(levels(test.data$Block))
              blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set

              test.data$blockeff =  blkeff

              for(g in as.vector(test.data$geno)){

                test.data[test.data$geno == g,"genoeff"] = genoUef[genoUef$geno == g, "V1"]
              }

              test.data$toteff = test.data[,"blockeff"]  + test.data[,"genoeff"]

              r = cbind(method, env, Trait, predictability = round(cor(test.data[,Trait],test.data[,"toteff"], use = "pairwise.complete.obs"),3))
              colnames(r)[4] <- c("predictability")
              method_outputs = rbind(method_outputs,r)

            }else if(method == "Block+Marker+Spatial"){
              # Design matrix
              #####################
              # make the dsisgn matrix for the blk_rwo and Blk-col
              # Random factor matrix
              dim(test.data)
              idblock <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
              Z.Blk.test <- model.matrix(~idblock - 1)

              rownames(Z.Blk.test) <- test.data$plot


              eval(parse(text = paste("ans<- mmer(",Trait,"~1,
                                     random=~ Block + vsr(geno, Gu = A1) +
                                       spl2Da(row,col),
                                     rcov=~vsr(units),
                                     data= train.data)")))

              ### Estimating effects
              befall =  as.matrix(ans$U$Block[[Trait]])
              len_b =as.numeric(levels(test.data$Block))
              blockeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set

              # make a plot to observe the spatial effects found by the spl2D()
              W <- with(test.data,spl2Da(row,col)) # 2D spline incidence matrix
              test.data$spatial <- W$Z$`A:all`%*%ans$U$`A:all`[[Trait]] # 2D spline BLUPs

              genoUef = as.matrix(ans$U$`u:geno`[[Trait]])
              genoUef = as.matrix(genoUef[order(rownames(genoUef)),])
              rownames(genoUef) %in% rownames(A1)
              id = which(rownames(genoUef) %in% unique(test.data$geno))
              genoUtest = as.matrix(genoUef[id,])

              rownames(genoUtest) %in% unique(test.data$geno)
              dim(genoUtest)
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
              efftest = blockeff + test.data$spatial + test.data$genoUtest
              test.data$toteff = efftest


              r = cbind(method,env,Trait,
                        predictability = round(cor(test.data[,paste(Trait)],
                                                   test.data[,"toteff"],
                                                   use = "pairwise.complete.obs"),3))
              colnames(r)[4] = "predictability"
              method_outputs = rbind(method_outputs,r)


            }
          }
        }
      }
    }

dim(method_outputs)
method_outputs[method_outputs$Trait ==  "Stripe.rust.plant.response...0.9.Mc.Neal.scale.CO_321.0001394",]
dim(method_outputs)
write.csv(x = dfcomb, file = "/Users/tas286/Documents/GitHub/Genomic selection wheat/spatialanaysis_final/sumpred_YLDQt_val_2014_Block.csv")
##########################
# ploting the graph

##########################
head(method_outputs)
Traits
unique(method_outputs$Trait)
acbac <- method_outputs[method_outputs$Trait %in% c("Stripe.rust.plant.response...0.9.Mc.Neal.scale.CO_321.0001394",
                        "Stripe.rust.severity.....CO_321.0001396"),]
acbacBlock <- acbac[acbac$method == "Block",]
dim(acbacBlock)
head(acbacBlock)
colnames(acbacBlock)[4] = "Block"
unique(acbac$method)
acbacBlsp <- acbac[acbac$method == "Block+Spatial" ,]
dim(acbacBlsp)
colnames(acbacBlsp)[4] = "Block+Spatial"
acbacBlMK <- acbac[acbac$method == "Block+Marker"  ,]
colnames(acbacBlMK)[4] = "Block+Marker"
acbacBlMKSp <- acbac[acbac$method == "Block+Marker+Spatial"  ,]
colnames(acbacBlMKSp)[4] <- "Block+Marker+Spatial"
acbac_comb <- data.frame(acbacBlock, acbacBlsp[,4], acbacBlMK[,4], acbacBlMKSp[,4])
colnames(acbac_comb)[5:7] <- c("Block_Spatial","Block_Marker","Block_Marker_Spatial")
str(acbac_comb)
acbac_comb$Block = as.numeric(acbac_comb$Block)
acbac_comb$Block_Spatial = as.numeric(acbac_comb$Block_Spatial)
acbac_comb$Block_Marker = as.numeric(acbac_comb$Block_Marker)
acbac_comb$Block_Marker_Spatial = as.numeric(acbac_comb$Block_Marker_Spatial)
p1 <- ggplot(data = acbac_comb, mapping = aes(x = Block, y = Block_Spatial, fill = env)) +
  geom_point(aes(colour = env))  + facet_grid(~Trait) +
  xlab("Block") + ylab("Block + Spatial")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw()
  # geom_point(data=Mean_ac_method,
  #            col="blue",
  #            size = 4,
  #            shape = 19,
  #            fill = "red")
p2 <- ggplot(data = acbac_comb, mapping = aes(x = Block, y = Block_Marker, fill = env)) +
  geom_point(aes(colour = env),)  +
  xlab("Block") + ylab("Block + Marker") + labs(title = "Bacterial Leaf Streak") + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw()
p3 <- ggplot(data = acbac_comb, mapping = aes(x = Block_Marker, y = Block_Marker_Spatial, fill = env)) +
  geom_point(aes(colour = env),)  + facet_grid(~Trait)+
  xlab("Block + Marker") + ylab("Block + Marker + Spatial") + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw()

library(ggpubr)
ggarrange(p1,p3, heights =8,legend = "bottom",common.legend = TRUE ,
          ncol = 1, nrow = 2, labels = c("A","B"),hjust = -2, vjust = 1) +
  theme(axis.title.y  = element_text(size = 1))






## stacked box plot
str(method_outputs)

dfcomb = method_outputs
dfcomb$predictability = as.numeric(dfcomb$predictability)
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




######################
# scatter plot for the differnt methods with and without spatial
library(readxl)
crosval_dat = read_excel("~/Documents/GitHub/Genomic_selection_wheat/spatialanaysis_final/The_corss_validation_saptial_final_with_Block.xlsx", sheet = 1)
head(crosval_dat)
crosval_Block = crosval_dat[crosval_dat$method == "Block",]
crosval_BlockSp = crosval_dat[crosval_dat$method == "Block+Spatial",]
crosval_BlockMark = crosval_dat[crosval_dat$method == "Block+Marker",]
crosval_BlockMarkSpat = crosval_dat[crosval_dat$method == "Block+Marker+Spatial",]

dim(crosval_BlockMark)

dim(crosval_BlockMarkSpat)


#########################
# read crossvalidation results -

paths = "~/Documents/GitHub/Genomic_selection_wheat/spatialanaysis_final"
list_files = c("sumpred_ABBmid_2014_Blacksburg_warsaw_Block.csv",
               "sumpred_ABBMid_2014_Lexington_woodford_Block.csv",
               "sumpred_CSR_Val_2015_Mead_Block.csv",
               "sumpred_HWWpanel_2012_Mead_Block.csv",
               "sumpred_HWWpanel_2012_Tipton_Block.csv",
               "sumpred_HWWpanel_2013_Tipton_Block.csv",
               "sumpred_YLDQt_Val_2014_Block.csv")
library(tibble)
cross_vl_result = tibble()
for(f in list_files){
  res1 = read.csv(file = paste(c(paths,f),collapse = "/"))
  cross_vl_result = rbind(cross_vl_result, res1)
}
unique(cross_vl_result$env)
cross_vl_result_Yield = cross_vl_result[!cross_vl_result$Trait == "Grain_Yield",  ]

cross_vl_result_complete = cross_vl_result[!cross_vl_result$env %in% c("HWWpanel_2012_Mead", "HWWpanel_2012_Tipton","HWWpanel_2013_Tipton" ,"YldQtl-Val_2014_ClayCenter", "YldQtl-Val_2014_Mead"),  ]
cross_vl_result_complete_A = cross_vl_result[cross_vl_result$env %in% c("HWWpanel_2012_Mead", "HWWpanel_2012_Tipton","HWWpanel_2013_Tipton" ,"YldQtl-Val_2014_ClayCenter", "YldQtl-Val_2014_Mead"),  ]
unique(cross_vl_result_complete_A$Trait)
table(cross_vl_result_complete$env,cross_vl_result_complete$Trait)
cross_vl_block = cross_vl_result_complete[cross_vl_result_complete$method == "Block",]
cross_vl_block = cross_vl_block[order(cross_vl_block$Trait),]
dim(cross_vl_block)
cross_vl_blockSp = cross_vl_result_complete[cross_vl_result_complete$method %in% "Block+Spatial",]
cross_vl_blockSp = cross_vl_blockSp[order(cross_vl_blockSp$Trait), ]
view(cross_vl_blockSp)

cross_vl_blockSp_sel = tibble()
for(trait in unique(cross_vl_blockSp$Trait)){
  for(env in unique(cross_vl_blockSp$env)){
    dt = cross_vl_blockSp[cross_vl_blockSp$env == env & cross_vl_blockSp$Trait == trait,]
    dt1 = dt[1:25,]
    cross_vl_blockSp_sel = rbind(cross_vl_blockSp_sel, dt1)
  }
}

dim(cross_vl_blockSp_sel)
cross_vl_blockSp_sel = cross_vl_blockSp_sel[order(cross_vl_blockSp_sel$Trait), ]
cross_vl_blockMark = cross_vl_result_complete[cross_vl_result_complete$method == "Block+Marker",]
cross_vl_blockMark = cross_vl_blockMark[order(cross_vl_blockMark$Trait),]
cross_vl_blockMarkSp = cross_vl_result_complete[cross_vl_result_complete$method == "Block+Marker+Spatial",]
cross_vl_blockMarkSp= cross_vl_blockMarkSp[order(cross_vl_blockMarkSp$Trait), ]
Cross_val_allMethod_wide = cbind(cross_vl_block,cross_vl_blockSp_sel,cross_vl_blockMark,cross_vl_blockMarkSp)
colnames(Cross_val_allMethod_wide)

####
# Trials with incomplet trait
cross_vl_result_complete_A = cross_vl_result[cross_vl_result$env %in% c("HWWpanel_2012_Mead", "HWWpanel_2012_Tipton","HWWpanel_2013_Tipton" ,"YldQtl-Val_2014_ClayCenter", "YldQtl-Val_2014_Mead"),  ]
unique(cross_vl_result_complete_A$Trait)
table(cross_vl_result_complete_A$env,cross_vl_result_complete_A$Trait)

cross_vl_Grain_Yield = cross_vl_result_complete_A[cross_vl_result_complete_A$Trait == "Grain_Yield",]
dim(cross_vl_Grain_Yield)
cross_vl_Grain_Yield_block = cross_vl_Grain_Yield[cross_vl_Grain_Yield$method == "Block",]
dim(cross_vl_Grain_Yield_block)
cross_vl_Grain_Yield_block = cross_vl_Grain_Yield_block[order(cross_vl_Grain_Yield_block$env), ]
table(cross_vl_Grain_Yield_block$env,cross_vl_Grain_Yield_block$Trait)
cross_vl_Grain_Yield_blocksp = cross_vl_Grain_Yield[cross_vl_Grain_Yield$method == "Block+Spatial",]
dim(cross_vl_Grain_Yield_blocksp)
cross_vl_Grain_Yield_blocksp = cross_vl_Grain_Yield_blocksp[order(cross_vl_Grain_Yield_blocksp$env),]
table(cross_vl_Grain_Yield_blocksp$env,cross_vl_Grain_Yield_blocksp$Trait)

cross_vl_Grain_Yield_blocksp_sel = tibble()
for(i in unique(cross_vl_Grain_Yield_blocksp$env)){
  for(j in unique(cross_vl_Grain_Yield_blocksp$method)){
    dg = cross_vl_Grain_Yield_blocksp[cross_vl_Grain_Yield_blocksp$env == i & cross_vl_Grain_Yield_blocksp$method == j, ]
    dg1 = dg[1:25,]
    cross_vl_Grain_Yield_blocksp_sel = rbind(cross_vl_Grain_Yield_blocksp_sel, dg1)
  }
}



cross_vl_Grain_Yield_blocMark = cross_vl_Grain_Yield[cross_vl_Grain_Yield$method == "Block+Marker",]
dim(cross_vl_Grain_Yield_blocMark)
cross_vl_Grain_Yield_blocMark = cross_vl_Grain_Yield_blocMark[order(cross_vl_Grain_Yield_blocMark$env),]
cross_vl_Grain_Yield_blocMarkSpt = cross_vl_Grain_Yield[cross_vl_Grain_Yield$method == "Block+Marker+Spatial",]
dim(cross_vl_Grain_Yield_blocMarkSpt)
Corsval_GY_Combined = cbind(cross_vl_Grain_Yield_block,cross_vl_Grain_Yield_blocksp_sel,cross_vl_Grain_Yield_blocMark,cross_vl_Grain_Yield_blocMarkSpt)
######
cross_vl_Plant_height = cross_vl_result_complete_A[cross_vl_result_complete_A$Trait %in% "Plant_height",]
str(cross_vl_Plant_height)
unique(cross_vl_Plant_height$method)
cross_vl_Plant_height_blk = cross_vl_Plant_height[cross_vl_Plant_height$method == "Block", ]
table(cross_vl_Plant_height_blk$env, cross_vl_Plant_height_blk$Trait)
cross_vl_Plant_height_blk = cross_vl_Plant_height_blk[order(cross_vl_Plant_height_blk$env), ]
cross_vl_Plant_height_blksp = cross_vl_Plant_height[cross_vl_Plant_height$method == "Block+Spatial", ]
table(cross_vl_Plant_height_blksp$env, cross_vl_Plant_height_blksp$Trait)
cross_vl_Plant_height_blksp = cross_vl_Plant_height_blksp[order(cross_vl_Plant_height_blksp$env),]
cross_vl_Plant_height_blksp_sel = tibble()
for(i in unique(cross_vl_Plant_height_blksp$env)){
  for(j in unique(cross_vl_Plant_height_blksp$method)){
    dk = cross_vl_Plant_height_blksp[cross_vl_Plant_height_blksp$env == i & cross_vl_Plant_height_blksp$method == j,]
    dk1 = dk[1:25,]
    cross_vl_Plant_height_blksp_sel = rbind(cross_vl_Plant_height_blksp_sel,dk1)
  }
}


cross_vl_Plant_height_blkMark = cross_vl_Plant_height[cross_vl_Plant_height$method == "Block+Marker", ]
table(cross_vl_Plant_height_blkMark$env, cross_vl_Plant_height_blkMark$Trait)
cross_vl_Plant_height_blkMark = cross_vl_Plant_height_blkMark[order(cross_vl_Plant_height_blkMark$env), ]
cross_vl_Plant_height_blkMarksp = cross_vl_Plant_height[cross_vl_Plant_height$method == "Block+Marker+Spatial", ]

table(cross_vl_Plant_height_blkMarksp$env, cross_vl_Plant_height_blkMarksp$Trait)
cross_vl_Plant_height_blkMarksp = cross_vl_Plant_height_blkMarksp[order(cross_vl_Plant_height_blkMarksp$env),]

cross_vl_Plant_height_combined = cbind(cross_vl_Plant_height_blk,cross_vl_Plant_height_blksp_sel,
                                       cross_vl_Plant_height_blkMark,cross_vl_Plant_height_blkMarksp)

cross_val_PH_GY = rbind(cross_vl_Plant_height_combined,Corsval_GY_Combined)
dim(cross_val_PH_GY)
dim(Cross_val_allMethod_wide)
cross_val_all_method_all_trait = rbind(cross_val_PH_GY,Cross_val_allMethod_wide )
cross_val_all_method_all_trait_sel = cross_val_all_method_all_trait[,c(1,2,3,7,11,15)]
colnames(cross_val_all_method_all_trait_sel) <- c("env", "Trait", "Acc_Block","Acc_Block_Spatial",
                                                  "Acc_Block_Marker","Acc_Block_Marker_Spatial" )
####
head(cross_val_all_method_all_trait_sel)

library(ggplot2)
library(dplyr)
Mean_ac_method = cross_val_all_method_all_trait_sel %>% group_by(Trait) %>%
  summarise_at(.vars = c("Acc_Block","Acc_Block_Spatial", "Acc_Block_Marker",
                         "Acc_Block_Marker_Spatial" ), .funs = mean, na.rm = T)

p1 = ggplot(data = cross_val_all_method_all_trait_sel, mapping = aes(x = Acc_Block, y = Acc_Block_Spatial, fill = env)) +
  geom_point(aes(colour = env))  + facet_wrap(~Trait)  +
  xlab("Block") + ylab("Block + Spatial")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw() +
  geom_point(data=Mean_ac_method,
             col="blue",
             size = 4,
             shape = 19,
             fill = "red")

p2 = ggplot(data = cross_val_all_method_all_trait_sel, mapping = aes(x = Acc_Block, y = Acc_Block_Marker)) +
  geom_point(aes(colour = env)) + facet_wrap(~Trait)  +
  xlab("Block") + ylab("Block + Marker")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw() +
  geom_point(data=Mean_ac_method,
             col="blue",
             size = 4,
             shape = 19,
             fill = "red")


p4 = ggplot(data = cross_val_all_method_all_trait_sel, mapping = aes(x = Acc_Block_Spatial, y = Acc_Block_Marker)) +
  geom_point(aes(colour = env)) + facet_wrap(~Trait) +
  xlab("Block + Spatial") + ylab("Block + Marker")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw() +
  geom_point(data=Mean_ac_method,
             col="blue",
             size = 4,
             shape = 19,
             fill = "red")


p3 = ggplot(data = cross_val_all_method_all_trait_sel, mapping = aes(x = Acc_Block, y = Acc_Block_Marker_Spatial)) +
  geom_point(aes(colour = env)) + facet_wrap(~Trait)  +
  xlab("Block") + ylab("Block + Marker + Spatial")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw() +
  geom_point(data=Mean_ac_method,
             col="blue",
             size = 4,
             shape = 19,
             fill = "red")

p5 = ggplot(data = cross_val_all_method_all_trait_sel, mapping = aes(x = Acc_Block_Spatial, y = Acc_Block_Marker_Spatial)) +
  geom_point(aes(colour = env)) + facet_wrap(~Trait) +
  xlab("Block + Spatial") + ylab("Block + Marker + Spatial")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw() +
  geom_point(data=Mean_ac_method,
             col="blue",
             size = 4,
             shape = 19,
             fill = "red")

p6 = ggplot(data = cross_val_all_method_all_trait_sel, mapping = aes(x = Acc_Block_Marker, y = Acc_Block_Marker_Spatial)) +
  geom_point(aes(colour = env)) + facet_wrap(~Trait)  +
  xlab("Block + Marker") + ylab("Block + Marker + Spatial")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw() +
  geom_point(data=Mean_ac_method,
             col="blue",
             size = 4,
             shape = 19,
             fill = "red")

library(gridExtra)
library(ggpubr)
p1 + theme(plot.margin = margin(10,8,10,8, "in"))
ggarrange(p1,p4,p2, heights =8,legend = "bottom",common.legend = TRUE ,
          ncol = 1, nrow = 3, labels = c("A","B","C"),hjust = -2, vjust = 2) +
  theme(axis.title.y  = element_text(size = 1))
ggarrange(p3 + ylab("") ,p5,p6 +  ylab(""), legend = "bottom", common.legend = TRUE,
          ncol = 1, nrow = 3,
          labels = c("A","B","C"),hjust = -2.5, vjust = 1.5) + theme(axis.title.y  = element_text(size = 3))

legend()
########
# plotting the predictability of the block vs block+ spatial for 52 trials
prd1 = read.csv("~/Documents/GitHub/Genomic_selection_wheat/50Trial_data/crossval_block_spat_24_trials.csv")
head(prd1)
prd1 = prd1[,-1]
prd2 = read.csv("~/Documents/GitHub/Genomic_selection_wheat/50Trial_data/mehtod_outputs_35_trials.csv")
head(prd2)
acc_58Trial <- rbind(prd1,prd2)
prd3 = read.csv("~/Documents/GitHub/Genomic_selection_wheat/50Trial_data/Method_accuracy_58_trials.csv")
dim(prd3)

head(acc_58Trial)
unique(acc_58Trial$method)
acc58_blk = acc_58Trial[acc_58Trial$method == "Block" , ]
dim(acc58_blk)

acc58_blkSp = acc_58Trial[acc_58Trial$method == "Block+Spatial" , ]
dim(acc58_blkSp)
acc58_comb_wide = cbind(acc58_blk,acc58_blkSp)
head(acc58_comb_wide)
colnames(acc58_comb_wide)[1] ="Block"
colnames(acc58_comb_wide)[4] ="Acc_Block"
colnames(acc58_comb_wide)[5] ="Block_Spatial"
colnames(acc58_comb_wide)[8] ="Acc_Block_Spatial"
acc58_comb_wide_sel = acc58_comb_wide[,c(1,2,3,4,5,8)]
library(ggplot2)

ggplot(data = acc58_comb_wide_sel, mapping = aes(x = Acc_Block, y = Acc_Block_Spatial)) +
  geom_point() + facet_wrap(~Trait) +geom_point(colour = "blue") +
  xlab("Block") + ylab("Block + Spatial")  + ylim(c(-0.4,1)) + xlim(c(-0.4,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw()


####Plotting heritability

hert = read.csv("~/Documents/GitHub/Genomic_selection_wheat/50Trial_data/Variance_comp_58_trials.csv")
head(hert)
ggplot(hert, mapping = aes(x = Trait, y = h2std, fill = method)) +
  geom_boxplot()
hertBl= hert[hert$method == "Block",]
hertBlSp= hert[hert$method == "Block+Spatial",]
hertcomb = cbind(hertBl,hertBlSp)
head(hertcomb)
colnames(hertcomb)[7] = "H2_Block"
colnames(hertcomb)[18] = "H2_Block_Spatial"
hertcomb_sel = hertcomb[, c(1,2,3,7,18)]
head(hertcomb_sel)

pbs = ggplot(data = hertcomb_sel, mapping = aes(x = H2_Block, y = H2_Block_Spatial)) +
  geom_point() + facet_wrap(~Trait) +geom_point(colour = "blue") +
  xlab("Block") + ylab("Block + Spatial")  + ylim(c(0,1)) + xlim(c(0,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw()

hert12 = read.csv("~/Documents/GitHub/Genomic_selection_wheat/12trials_Marker-data/Hertiability_estimate_12_trials.csv")
head(hert12)
hert12_MSp = hert12[hert12$Method == "W_m_w_sp",]
hert12_M = hert12[hert12$Method == "W_m_wo_sp",]
hert12Comb_MSP = cbind(hert12_MSp,hert12_M)
head(hert12Comb_MSP)
colnames(hert12Comb_MSP)[4] = "H2_M_Spatial"
colnames(hert12Comb_MSP)[8] = "H2_Mark"
hert12Comb_MSP_sel = hert12Comb_MSP[,c(2,3,4,8)]

pms = ggplot(data = hert12Comb_MSP_sel, mapping = aes(x = H2_Mark, y = H2_M_Spatial)) +
  geom_point() + facet_wrap(~Trait) +geom_point(colour = "blue") +
  xlab("Marker") + ylab("Marker+ Spatial")  + ylim(c(0,1)) + xlim(c(0,1)) +
  geom_abline(intercept = 0, colour = "red", linewidth = 1.2) + theme_bw()
library(ggpubr)
ggarrange(pms, pbs,ncol = 1, nrow = 2,labels = c("A","B"),
          font.label = list(size = 12, face = "plain"),
          hjust = -4,vjust = 2)
