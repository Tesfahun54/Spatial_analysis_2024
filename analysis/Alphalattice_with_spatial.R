AlphaLattice_with_spatial <- function(Trait, loc_data){
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
return(output)
}
