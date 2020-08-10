# PRS_comparisons
Using AUC analysis to see how well the polygenic risk profiles of RBD and of various conditions can differentiate between RBD cases and controls.  

Polygenic risk profiles were defined as GWAS significant (and replicated) variants, except for DLB, were the p-value threshold was determined by permutation by Scholz et. al., and RBD, where the FDR threshold was applied. Risk scores were calulcated in an idiopathic RBD cohort (cases=1037, controls=909) using plink.  

Genetic profiles to test were determined based on previous associations to PD and/or genetic correlation with RBD. 

## Testing RBD PRS performance in RBD
```R
require(data.table)
require(pROC)

RBD_FDR <- fread("rbd_meta_fdr_adjusted.profile", header = T)
RBD_FDR <- subset(RBD_FDR, PHENO != -9)
RBD_FDR$PHENOT <-ifelse(RBD_FDR$PHENO==1, "CONTROL", ifelse(RBD_FDR$PHENO==2,"RBD", NA))

rocAuc <- roc(RBD_FDR$PHENO, RBD_FDR$SCORE)
auc(rocAuc) # 0.662
ci(rocAuc, of="auc") # 0.637-0.685 (DeLong)
coords(rocAuc, "best") # threshold: 0.00266 # specificity: 0.60 # sensitivity: 0.66
thresh <- coords(rocAuc, "best")[1]

png("TEST_AUC_RBD-FDR_cont.png", width = 5.5, height = 4, units = "in", res = 300)

rocobj <- plot.roc(RBD_FDR$PHENO, RBD_FDR$SCORE,  main="ROC curve RBD vs controls",
                   percent=FALSE,  ci=TRUE, print.auc=TRUE, col = "darkred")
ciobj <- ci.se(rocobj, specificities=seq(0, 1, 0.05))
plot(ciobj, type="shape", col="gray90")

dev.off()
```
![RBD FDR AUC](TEST_AUC_RBD-FDR_cont.png)


## Testing genetically correlated conditions 
```R

SCZ <- fread("scz_test.profile", header = T) 
BIP <- fread("bipolar_test.profile", header = T)
INS <- fread("insomnia_test.profile", header = T)
NARC <- fread("narc_test.profile", header = T)
MS <- fread("MS_test.profile", header = T)
AD <- fread("AD_test.profile", header = T)
PD <- fread("PD_test.profile", header = T)
DIAB2 <- fread("diab2_test.profile", header = T)

DLB <- fread("dlb_sonja_hg19.profile", header = T)
DLB <- subset(DLB, PHENO != -9)
DLB_noAPOE <- fread("dlb_sonja_hg19_noAPOE.profile", header = T)
DLB_noAPOE <- subset(DLB_noAPOE, PHENO != -9)
````

## Create ROC curve
```R
png("TEST_AUC_RBD_cont.png", width = 5.5, height = 4, units = "in", res = 300)

rocobj <- plot.roc(SCZ$PHENOT, SCZ$SCORE,  main="ROC curve RBD vs controls: Correlated Conditions",
                   percent=FALSE,  ci=TRUE, print.auc=F, col = "darkred")
rocobj <- plot(roc(INS$PHENOT, INS$SCORE), print.auc = F, 
     col = "green", add = TRUE)
rocobj <- plot(roc(BIP$PHENOT, BIP$SCORE), print.auc = F, 
               col = "blue", add = TRUE)
rocobj <- plot(roc(DIAB2$PHENOT, DIAB2$SCORE), print.auc = F, 
                col = "darkorange", add = TRUE)
rocobj <- plot(roc(NARC$PHENOT, NARC$SCORE), print.auc = F, 
               col = "darkgreen", add = TRUE)
rocobj <- plot.roc(RBD$PHENOT, RBD$SCORE,  main="ROC curve RBD vs controls",
                   percent=FALSE,  ci=TRUE, print.auc=TRUE, col = "darkblue", add = TRUE)

rocobj <- plot.roc(AD$PHENOT, AD$SCORE,  main="ROC curve RBD vs controls",
                    percent=FALSE,  ci=TRUE, print.auc=F, col = "red", add = TRUE)
rocobj <- plot(roc(PD$PHENOT, PD$SCORE), print.auc = F, 
               col = "darkblue", add = TRUE)
rocobj <- plot(roc(MS$PHENOT, MS$SCORE), print.auc = F, 
                       col = "darkorange", add = TRUE)
rocobj <- plot(roc(DIAB2$PHENOT, DIAB2$SCORE), print.auc = F, 
               col = "brown", add = TRUE)

dev.off()
```
![ROC Results](TEST_AUC_RBD_cont.png)
