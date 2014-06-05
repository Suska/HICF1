HICF1 - factor modelling - TREES - V2
========================================================

Dr. Susanne Weller 
27/05/2014

Import merged dataframe - v2


```r
setwd("~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v3")
load("genclinv3.Rda")
library(plyr)

# reorder collumns for convenience
genclinv3 <- genclinv3[c(1, 4, 8, 2, 6, 5, 50, 3, 7, 9:11, 13, 15:50, 53:57)]
genclinv3$X11q_monoallelic_mut <- NULL
genclinv3$Binet.x.1 <- NULL
genclinv3 <- rename(genclinv3, c(X11q_monoallelic_del.y = "X11q_mono_del", X11q_biallelic.y = "X11q_biallelic"))
names(genclinv3)[7] <- "Binet"

# convert integers to factors

genclinv3$vhmut <- as.factor(genclinv3$vhmut)
levels(genclinv3$MRD) <- c("-MRD", "+MRD")
genclinv3$cd38 <- as.factor(genclinv3$cd38)
genclinv3[, 13:47] <- lapply(genclinv3[, 13:47], as.factor)

# combine genetic factors
genclinv3$SAMHD1 <- factor(ifelse(genclinv3$SAMHD1_1mutationONLY == "1" | genclinv3$SAMHD1_biallelic_all == 
    "1", "1", "0"))
genclinv3$TP53 <- factor(ifelse(genclinv3$TP53_mutationONLY == "1" | genclinv3$TP53biallelic == 
    "1", "1", "0"))
genclinv3 <- genclinv3[c(1:47, 53:54, 48:52)]
```


SIMPLE CLASSIFICATION TREES - with blood data added
===================================================


```r
treegenclin <- genclinv3[c(6:54)]

library(rpart)
library(rpart.plot)
source("lastfunction.R")

# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v3/figures_trees_03062014/01_MissClassErr_alldata.pdf')
par(mfrow = c(1, 1))
for (cp in c(0.01, 0.009, 0.008, 0.006, 0.004, 0.002, 0.001)) {
    source("lastfunction.R")
    nobuckets <- NULL
    MissClassErr <- NULL
    CorrectClass <- NULL
    rootnoderr <- NULL
    relErr <- NULL
    for (i in 1:15) {
        # built the tree
        rp.treegenclin <- rpart(treegenclin$MRD ~ ., method = "class", data = treegenclin, 
            control = rpart.control(minbucket = i, xval = 10, cp = cp))
        # get the number of buckets
        nobuckets[i] <- i
        # get the root node error
        rootnoderr[i] <- rp.treegenclin$parms$prior[2]
        # Get the rel error for each tree
        rpsum <- printcp(rp.treegenclin)
        relErr[i] <- last(rpsum)
        # Calculate missclassification error
        MissClassErr[i] <- rootnoderr[i] * relErr[i]
        # Calculate correctly classified patients
        CorrectClass[i] <- 209 - MissClassErr[i] * 209
    }
    MissClassErr
    if (cp == 0.01) {
        plot(nobuckets, MissClassErr, pch = 0, type = "o", ylim = c(0, 0.4), 
            main = "All data without age, treatment, freqvhmut", xlab = "No. of patients in end group\n(more patients=less complex model)", 
            ylab = "Missclassification error=fraction of patients that will be put in wrong group ")
        legend(12, 0.2, c("0.01", "0.009", "0.008", "0.006", "0.004", "0.002", 
            "0.001"), cex = 0.8, col = c("black", "green", "blue", "orange", 
            "purple", "darkslategray4", "red"), pch = 0:6, title = "complexity:\n lower = more complex")
    }
    if (cp == 0.009) {
        lines(x = nobuckets + 0.05, MissClassErr, type = "o", pch = 1, col = "green")
    }
    if (cp == 0.008) {
        lines(x = nobuckets + 0.1, MissClassErr, type = "o", pch = 2, col = "blue")
    }
    if (cp == 0.006) {
        lines(x = nobuckets + 0.15, MissClassErr, type = "o", pch = 3, col = "orange")
    }
    if (cp == 0.004) {
        lines(x = nobuckets + 0.2, MissClassErr, type = "o", pch = 4, col = "purple")
    }
    if (cp == 0.002) {
        lines(x = nobuckets + 0.25, MissClassErr, type = "o", pch = 5, col = "darkslategray4")
    }
    if (cp == 0.001) {
        lines(x = nobuckets + 0.3, MissClassErr, type = "o", pch = 6, col = "red")
    }
}
```

```
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   0.95 0.071
## 3 0.049      4      0.61   0.87 0.070
## 4 0.039      5      0.56   0.84 0.070
## 5 0.027      7      0.48   0.78 0.069
## 6 0.020     12      0.34   0.75 0.068
## 7 0.015     19      0.21   0.72 0.068
## 8 0.010     21      0.18   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.06 0.071
## 2 0.065      1      0.80   0.87 0.070
## 3 0.049      4      0.61   0.79 0.069
## 4 0.039      5      0.56   0.75 0.068
## 5 0.027      7      0.48   0.71 0.067
## 6 0.020     12      0.34   0.63 0.065
## 7 0.013     19      0.21   0.60 0.064
## 8 0.010     22      0.17   0.60 0.064
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           plat           SAMHD1        
##  [9] TP53           WBC            X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.065      1      0.80   0.83 0.070
## 3 0.049      4      0.61   0.85 0.070
## 4 0.039      5      0.56   0.82 0.069
## 5 0.027      7      0.48   0.68 0.067
## 6 0.020     12      0.34   0.67 0.066
## 7 0.013     14      0.30   0.66 0.066
## 8 0.010     17      0.26   0.65 0.066
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] trisomy12      vhmut          WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.92 0.071
## 3 0.049      4      0.61   0.77 0.069
## 4 0.039      5      0.56   0.75 0.068
## 5 0.025      7      0.48   0.75 0.068
## 6 0.020      9      0.43   0.81 0.069
## 7 0.013     14      0.33   0.82 0.069
## 8 0.010     17      0.29   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] plat           SAMHD1         TP53           WBC           
## [9] X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   0.94 0.071
## 3 0.049      4      0.61   0.82 0.069
## 4 0.039      7      0.46   0.79 0.069
## 5 0.029      8      0.42   0.77 0.069
## 6 0.026     10      0.36   0.77 0.069
## 7 0.010     13      0.28   0.82 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.11 0.071
## 2 0.065      1      0.80   0.88 0.070
## 3 0.049      4      0.61   0.82 0.069
## 4 0.029      7      0.46   0.75 0.068
## 5 0.026      8      0.43   0.73 0.068
## 6 0.010     11      0.35   0.73 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.15 0.070
## 2 0.065      1      0.80   0.93 0.071
## 3 0.049      4      0.61   0.87 0.070
## 4 0.029      7      0.46   0.70 0.067
## 5 0.015      8      0.43   0.72 0.068
## 6 0.010     10      0.40   0.73 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.065      1      0.80   0.98 0.071
## 3 0.049      4      0.61   0.90 0.070
## 4 0.039      6      0.51   0.89 0.070
## 5 0.029      7      0.47   0.79 0.069
## 6 0.010      8      0.44   0.81 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.92 0.071
## 3 0.049      4      0.61   0.89 0.070
## 4 0.025      6      0.51   0.86 0.070
## 5 0.010      8      0.46   0.87 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   0.89 0.070
## 3 0.049      4      0.61   0.78 0.069
## 4 0.020      6      0.51   0.80 0.069
## 5 0.010      8      0.47   0.81 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.90 0.070
## 3 0.049      4      0.61   0.88 0.070
## 4 0.020      6      0.51   0.78 0.069
## 5 0.010      8      0.47   0.79 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   0.83 0.070
## 3 0.049      4      0.61   0.73 0.068
## 4 0.020      6      0.51   0.70 0.067
## 5 0.010      8      0.47   0.72 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.94 0.071
## 3 0.049      4      0.61   0.88 0.070
## 4 0.020      6      0.51   0.81 0.069
## 5 0.010      8      0.47   0.81 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.18 0.070
## 2 0.083      1      0.80   0.87 0.070
## 3 0.039      3      0.64   0.79 0.069
## 4 0.020      6      0.52   0.85 0.070
## 5 0.010      7      0.50   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.083      1      0.80   0.91 0.070
## 3 0.039      3      0.64   0.84 0.070
## 4 0.020      6      0.52   0.90 0.070
## 5 0.010      7      0.50   0.88 0.070
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-21.png) 

```
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     Binet          BIRC3_mutation
##  [5] cd38           CNAs           del8p          haem          
##  [9] MED12mutation  plat           SAMHD1         SF3B1_mutation
## [13] TP53           TP53biallelic  WBC            X11q_biallelic
## [17] X11q_mono_del  X13q_hom      
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0     1.000   1.11 0.071
## 2 0.0654      1     0.804   0.86 0.070
## 3 0.0490      4     0.608   0.81 0.069
## 4 0.0392      5     0.559   0.80 0.069
## 5 0.0275      7     0.480   0.81 0.069
## 6 0.0196     12     0.343   0.80 0.069
## 7 0.0147     19     0.206   0.80 0.069
## 8 0.0098     21     0.176   0.78 0.069
## 9 0.0090     29     0.098   0.78 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   0.86 0.070
## 3 0.0490      4      0.61   0.79 0.069
## 4 0.0392      5      0.56   0.79 0.069
## 5 0.0275      7      0.48   0.77 0.069
## 6 0.0196     12      0.34   0.68 0.067
## 7 0.0131     19      0.21   0.69 0.067
## 8 0.0098     22      0.17   0.70 0.067
## 9 0.0090     26      0.13   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           plat           SAMHD1        
##  [9] SF3B1_mutation TP53           WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.11 0.071
## 2 0.0654      1      0.80   0.89 0.070
## 3 0.0490      4      0.61   0.75 0.068
## 4 0.0392      5      0.56   0.76 0.069
## 5 0.0275      7      0.48   0.76 0.069
## 6 0.0196     12      0.34   0.70 0.067
## 7 0.0131     14      0.30   0.76 0.069
## 8 0.0098     17      0.26   0.75 0.068
## 9 0.0090     20      0.24   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] trisomy12      vhmut          WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.18 0.070
## 2 0.065      1      0.80   0.96 0.071
## 3 0.049      4      0.61   0.85 0.070
## 4 0.039      5      0.56   0.70 0.067
## 5 0.025      7      0.48   0.67 0.066
## 6 0.020      9      0.43   0.67 0.066
## 7 0.013     14      0.33   0.67 0.066
## 8 0.009     17      0.29   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] plat           SAMHD1         TP53           WBC           
## [9] X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.065      1      0.80   0.96 0.071
## 3 0.049      4      0.61   0.81 0.069
## 4 0.039      7      0.46   0.78 0.069
## 5 0.029      8      0.42   0.77 0.069
## 6 0.026     10      0.36   0.74 0.068
## 7 0.009     13      0.28   0.69 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   0.93 0.071
## 3 0.0490      4      0.61   0.86 0.070
## 4 0.0294      7      0.46   0.73 0.068
## 5 0.0261      8      0.43   0.72 0.068
## 6 0.0098     11      0.35   0.68 0.067
## 7 0.0090     15      0.31   0.67 0.066
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.84 0.070
## 3 0.0490      4      0.61   0.82 0.069
## 4 0.0294      7      0.46   0.71 0.067
## 5 0.0147      8      0.43   0.72 0.068
## 6 0.0098     10      0.40   0.72 0.068
## 7 0.0090     14      0.36   0.71 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.95 0.071
## 3 0.0490      4      0.61   0.74 0.068
## 4 0.0392      6      0.51   0.76 0.069
## 5 0.0294      7      0.47   0.74 0.068
## 6 0.0098      8      0.44   0.77 0.069
## 7 0.0090     10      0.42   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.92 0.071
## 3 0.049      4      0.61   0.91 0.070
## 4 0.025      6      0.51   0.83 0.070
## 5 0.009      8      0.46   0.83 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.85 0.070
## 3 0.049      4      0.61   0.78 0.069
## 4 0.020      6      0.51   0.78 0.069
## 5 0.009      8      0.47   0.77 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.92 0.071
## 3 0.049      4      0.61   0.84 0.070
## 4 0.020      6      0.51   0.76 0.069
## 5 0.009      8      0.47   0.79 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.90 0.070
## 3 0.049      4      0.61   0.82 0.069
## 4 0.020      6      0.51   0.88 0.070
## 5 0.009      8      0.47   0.82 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.08 0.071
## 2 0.065      1      0.80   1.00 0.071
## 3 0.049      4      0.61   0.77 0.069
## 4 0.020      6      0.51   0.83 0.070
## 5 0.009      8      0.47   0.82 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.083      1      0.80   0.81 0.069
## 3 0.039      3      0.64   0.72 0.068
## 4 0.020      6      0.52   0.75 0.068
## 5 0.009      7      0.50   0.78 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.083      1      0.80   0.83 0.070
## 3 0.039      3      0.64   0.73 0.068
## 4 0.020      6      0.52   0.75 0.068
## 5 0.009      7      0.50   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     Binet          BIRC3_mutation
##  [5] cd38           CNAs           del8p          haem          
##  [9] MED12mutation  plat           SAMHD1         SF3B1_mutation
## [13] TP53           TP53biallelic  WBC            X11q_biallelic
## [17] X11q_mono_del  X13q_hom      
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0     1.000   1.07 0.071
## 2 0.0654      1     0.804   0.86 0.070
## 3 0.0490      4     0.608   0.82 0.069
## 4 0.0392      5     0.559   0.75 0.068
## 5 0.0275      7     0.480   0.69 0.067
## 6 0.0196     12     0.343   0.70 0.067
## 7 0.0147     19     0.206   0.72 0.068
## 8 0.0098     21     0.176   0.71 0.067
## 9 0.0080     29     0.098   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   1.00 0.071
## 3 0.0490      4      0.61   0.86 0.070
## 4 0.0392      5      0.56   0.86 0.070
## 5 0.0275      7      0.48   0.77 0.069
## 6 0.0196     12      0.34   0.71 0.067
## 7 0.0131     19      0.21   0.69 0.067
## 8 0.0098     22      0.17   0.69 0.067
## 9 0.0080     26      0.13   0.63 0.065
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           plat           SAMHD1        
##  [9] SF3B1_mutation TP53           WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   0.88 0.070
## 3 0.0490      4      0.61   0.75 0.068
## 4 0.0392      5      0.56   0.73 0.068
## 5 0.0275      7      0.48   0.70 0.067
## 6 0.0196     12      0.34   0.69 0.067
## 7 0.0131     14      0.30   0.71 0.067
## 8 0.0098     17      0.26   0.75 0.068
## 9 0.0080     20      0.24   0.77 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] trisomy12      vhmut          WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.98 0.071
## 3 0.049      4      0.61   0.92 0.071
## 4 0.039      5      0.56   0.86 0.070
## 5 0.025      7      0.48   0.69 0.067
## 6 0.020      9      0.43   0.72 0.068
## 7 0.013     14      0.33   0.70 0.067
## 8 0.008     17      0.29   0.74 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] plat           SAMHD1         TP53           WBC           
## [9] X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.065      1      0.80   0.95 0.071
## 3 0.049      4      0.61   0.94 0.071
## 4 0.039      7      0.46   0.83 0.070
## 5 0.029      8      0.42   0.82 0.069
## 6 0.026     10      0.36   0.82 0.069
## 7 0.008     13      0.28   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.91 0.070
## 3 0.0490      4      0.61   0.77 0.069
## 4 0.0294      7      0.46   0.65 0.066
## 5 0.0261      8      0.43   0.67 0.066
## 6 0.0098     11      0.35   0.71 0.067
## 7 0.0080     15      0.31   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.86 0.070
## 3 0.0490      4      0.61   0.77 0.069
## 4 0.0294      7      0.46   0.76 0.069
## 5 0.0147      8      0.43   0.78 0.069
## 6 0.0098     10      0.40   0.77 0.069
## 7 0.0080     14      0.36   0.80 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.94 0.071
## 3 0.0490      4      0.61   0.77 0.069
## 4 0.0392      6      0.51   0.72 0.068
## 5 0.0294      7      0.47   0.69 0.067
## 6 0.0098      8      0.44   0.74 0.068
## 7 0.0080     10      0.42   0.71 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.95 0.071
## 3 0.049      4      0.61   0.82 0.069
## 4 0.025      6      0.51   0.77 0.069
## 5 0.008      8      0.46   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.90 0.070
## 3 0.049      4      0.61   0.79 0.069
## 4 0.020      6      0.51   0.82 0.069
## 5 0.008      8      0.47   0.90 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.065      1      0.80   0.90 0.070
## 3 0.049      4      0.61   0.86 0.070
## 4 0.020      6      0.51   0.85 0.070
## 5 0.008      8      0.47   0.83 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.85 0.070
## 3 0.049      4      0.61   0.79 0.069
## 4 0.020      6      0.51   0.82 0.069
## 5 0.008      8      0.47   0.83 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.91 0.070
## 3 0.049      4      0.61   0.77 0.069
## 4 0.020      6      0.51   0.78 0.069
## 5 0.008      8      0.47   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.083      1      0.80   1.01 0.071
## 3 0.039      3      0.64   0.83 0.070
## 4 0.020      6      0.52   0.79 0.069
## 5 0.008      7      0.50   0.80 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.083      1      0.80   1.01 0.071
## 3 0.039      3      0.64   0.81 0.069
## 4 0.020      6      0.52   0.82 0.069
## 5 0.008      7      0.50   0.82 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     Binet          BIRC3_mutation
##  [5] cd38           CNAs           del8p          haem          
##  [9] MED12mutation  plat           SAMHD1         SF3B1_mutation
## [13] TP53           TP53biallelic  vhmut          WBC           
## [17] X11q_biallelic X11q_mono_del  X13q_hom      
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##        CP nsplit rel error xerror  xstd
## 1  0.1961      0     1.000   1.08 0.071
## 2  0.0654      1     0.804   0.99 0.071
## 3  0.0490      4     0.608   0.88 0.070
## 4  0.0392      5     0.559   0.81 0.069
## 5  0.0275      7     0.480   0.76 0.069
## 6  0.0196     12     0.343   0.73 0.068
## 7  0.0147     19     0.206   0.69 0.067
## 8  0.0098     21     0.176   0.69 0.067
## 9  0.0074     29     0.098   0.75 0.068
## 10 0.0065     33     0.069   0.75 0.068
## 11 0.0060     36     0.049   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.16 0.070
## 2 0.0654      1      0.80   0.94 0.071
## 3 0.0490      4      0.61   0.84 0.070
## 4 0.0392      5      0.56   0.76 0.069
## 5 0.0275      7      0.48   0.71 0.067
## 6 0.0196     12      0.34   0.66 0.066
## 7 0.0131     19      0.21   0.69 0.067
## 8 0.0098     22      0.17   0.71 0.067
## 9 0.0060     26      0.13   0.74 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           plat           SAMHD1        
##  [9] SF3B1_mutation TP53           WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   0.98 0.071
## 3 0.0490      4      0.61   0.73 0.068
## 4 0.0392      5      0.56   0.70 0.067
## 5 0.0275      7      0.48   0.75 0.068
## 6 0.0196     12      0.34   0.71 0.067
## 7 0.0131     14      0.30   0.74 0.068
## 8 0.0098     17      0.26   0.75 0.068
## 9 0.0060     20      0.24   0.81 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] trisomy12      vhmut          WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.87 0.070
## 3 0.049      4      0.61   0.80 0.069
## 4 0.039      5      0.56   0.72 0.068
## 5 0.025      7      0.48   0.67 0.066
## 6 0.020      9      0.43   0.66 0.066
## 7 0.013     14      0.33   0.71 0.067
## 8 0.006     17      0.29   0.69 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] plat           SAMHD1         TP53           WBC           
## [9] X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   0.93 0.071
## 3 0.049      4      0.61   0.77 0.069
## 4 0.039      7      0.46   0.79 0.069
## 5 0.029      8      0.42   0.71 0.067
## 6 0.026     10      0.36   0.67 0.066
## 7 0.006     13      0.28   0.69 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.87 0.070
## 3 0.0490      4      0.61   0.79 0.069
## 4 0.0294      7      0.46   0.74 0.068
## 5 0.0261      8      0.43   0.75 0.068
## 6 0.0098     11      0.35   0.79 0.069
## 7 0.0060     15      0.31   0.80 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.20 0.070
## 2 0.0654      1      0.80   0.98 0.071
## 3 0.0490      4      0.61   0.85 0.070
## 4 0.0294      7      0.46   0.78 0.069
## 5 0.0147      8      0.43   0.80 0.069
## 6 0.0098     10      0.40   0.82 0.069
## 7 0.0060     14      0.36   0.86 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho       ANC_Neutro       CNAs             haem            
## [5] independent_13q. SF3B1_mutation   TP53             WBC             
## [9] X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   0.92 0.071
## 3 0.0490      4      0.61   0.73 0.068
## 4 0.0392      6      0.51   0.77 0.069
## 5 0.0294      7      0.47   0.75 0.068
## 6 0.0098      8      0.44   0.75 0.068
## 7 0.0065     10      0.42   0.75 0.068
## 8 0.0060     13      0.40   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.10 0.071
## 2 0.065      1      0.80   0.97 0.071
## 3 0.049      4      0.61   0.87 0.070
## 4 0.025      6      0.51   0.83 0.070
## 5 0.006      8      0.46   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   0.82 0.069
## 3 0.049      4      0.61   0.78 0.069
## 4 0.020      6      0.51   0.75 0.068
## 5 0.006      8      0.47   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.08 0.071
## 2 0.065      1      0.80   0.94 0.071
## 3 0.049      4      0.61   0.79 0.069
## 4 0.020      6      0.51   0.75 0.068
## 5 0.006      8      0.47   0.77 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.15 0.070
## 2 0.065      1      0.80   0.84 0.070
## 3 0.049      4      0.61   0.93 0.071
## 4 0.020      6      0.51   0.90 0.070
## 5 0.006      8      0.47   0.83 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.90 0.070
## 3 0.049      4      0.61   0.71 0.067
## 4 0.020      6      0.51   0.70 0.067
## 5 0.006      8      0.47   0.67 0.066
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.13 0.071
## 2 0.083      1      0.80   0.98 0.071
## 3 0.039      3      0.64   0.85 0.070
## 4 0.020      6      0.52   0.88 0.070
## 5 0.006      7      0.50   0.86 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.16 0.070
## 2 0.083      1      0.80   0.86 0.070
## 3 0.039      3      0.64   0.75 0.068
## 4 0.020      6      0.52   0.86 0.070
## 5 0.006      7      0.50   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     Binet          BIRC3_mutation
##  [5] cd38           CNAs           del8p          gender        
##  [9] haem           MED12mutation  plat           SAMHD1        
## [13] SF3B1_mutation TP53           TP53biallelic  trisomy12     
## [17] vhmut          WBC            X11q_biallelic X11q_mono_del 
## [21] X13q_hom      
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##        CP nsplit rel error xerror  xstd
## 1  0.1961      0     1.000   1.06 0.071
## 2  0.0654      1     0.804   0.91 0.070
## 3  0.0490      4     0.608   0.85 0.070
## 4  0.0392      5     0.559   0.82 0.069
## 5  0.0275      7     0.480   0.82 0.069
## 6  0.0196     12     0.343   0.77 0.069
## 7  0.0147     19     0.206   0.79 0.069
## 8  0.0098     21     0.176   0.73 0.068
## 9  0.0074     29     0.098   0.83 0.070
## 10 0.0065     33     0.069   0.85 0.070
## 11 0.0049     36     0.049   0.85 0.070
## 12 0.0040     40     0.029   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.18 0.070
## 2 0.0654      1      0.80   0.95 0.071
## 3 0.0490      4      0.61   0.89 0.070
## 4 0.0392      5      0.56   0.82 0.069
## 5 0.0275      7      0.48   0.76 0.069
## 6 0.0196     12      0.34   0.74 0.068
## 7 0.0131     19      0.21   0.73 0.068
## 8 0.0098     22      0.17   0.72 0.068
## 9 0.0040     26      0.13   0.69 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           plat           SAMHD1        
##  [9] SF3B1_mutation TP53           WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##        CP nsplit rel error xerror  xstd
## 1  0.1961      0      1.00   1.00 0.071
## 2  0.0654      1      0.80   0.81 0.069
## 3  0.0490      4      0.61   0.80 0.069
## 4  0.0392      5      0.56   0.75 0.068
## 5  0.0275      7      0.48   0.67 0.066
## 6  0.0196     12      0.34   0.67 0.066
## 7  0.0131     14      0.30   0.65 0.066
## 8  0.0098     17      0.26   0.66 0.066
## 9  0.0049     20      0.24   0.71 0.067
## 10 0.0040     22      0.23   0.71 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] trisomy12      vhmut          WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   1.00 0.071
## 3 0.0490      4      0.61   0.83 0.070
## 4 0.0392      5      0.56   0.77 0.069
## 5 0.0245      7      0.48   0.73 0.068
## 6 0.0196      9      0.43   0.69 0.067
## 7 0.0131     14      0.33   0.70 0.067
## 8 0.0049     17      0.29   0.73 0.068
## 9 0.0040     19      0.28   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.06 0.071
## 2 0.0654      1      0.80   0.88 0.070
## 3 0.0490      4      0.61   0.91 0.070
## 4 0.0392      7      0.46   0.78 0.069
## 5 0.0294      8      0.42   0.70 0.067
## 6 0.0261     10      0.36   0.69 0.067
## 7 0.0049     13      0.28   0.67 0.066
## 8 0.0040     15      0.27   0.73 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.16 0.070
## 2 0.0654      1      0.80   1.00 0.071
## 3 0.0490      4      0.61   0.88 0.070
## 4 0.0294      7      0.46   0.67 0.066
## 5 0.0261      8      0.43   0.69 0.067
## 6 0.0098     11      0.35   0.58 0.064
## 7 0.0040     15      0.31   0.59 0.064
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.11 0.071
## 2 0.0654      1      0.80   0.84 0.070
## 3 0.0490      4      0.61   0.76 0.069
## 4 0.0294      7      0.46   0.69 0.067
## 5 0.0147      8      0.43   0.76 0.069
## 6 0.0098     10      0.40   0.80 0.069
## 7 0.0040     14      0.36   0.79 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho       ANC_Neutro       CNAs             haem            
## [5] independent_13q. SF3B1_mutation   TP53             WBC             
## [9] X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.86 0.070
## 3 0.0490      4      0.61   0.75 0.068
## 4 0.0392      6      0.51   0.77 0.069
## 5 0.0294      7      0.47   0.70 0.067
## 6 0.0098      8      0.44   0.73 0.068
## 7 0.0065     10      0.42   0.74 0.068
## 8 0.0040     13      0.40   0.74 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.88 0.070
## 3 0.049      4      0.61   0.82 0.069
## 4 0.025      6      0.51   0.76 0.069
## 5 0.004      8      0.46   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.065      1      0.80   1.03 0.071
## 3 0.049      4      0.61   0.87 0.070
## 4 0.020      6      0.51   0.76 0.069
## 5 0.004      8      0.47   0.74 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.11 0.071
## 2 0.065      1      0.80   0.98 0.071
## 3 0.049      4      0.61   0.93 0.071
## 4 0.020      6      0.51   0.84 0.070
## 5 0.004      8      0.47   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.95 0.071
## 3 0.049      4      0.61   0.83 0.070
## 4 0.020      6      0.51   0.87 0.070
## 5 0.004      8      0.47   0.80 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.94 0.071
## 3 0.049      4      0.61   0.87 0.070
## 4 0.020      6      0.51   0.86 0.070
## 5 0.004      8      0.47   0.82 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror xstd
## 1 0.196      0      1.00   1.14 0.07
## 2 0.083      1      0.80   0.91 0.07
## 3 0.039      3      0.64   0.83 0.07
## 4 0.020      6      0.52   0.89 0.07
## 5 0.004      7      0.50   0.84 0.07
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.11 0.071
## 2 0.083      1      0.80   0.97 0.071
## 3 0.039      3      0.64   0.79 0.069
## 4 0.020      6      0.52   0.84 0.070
## 5 0.004      7      0.50   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     Binet          BIRC3_mutation
##  [5] cd38           clones         CNAs           del4p         
##  [9] del8p          gender         haem           MED12mutation 
## [13] plat           SAMHD1         SF3B1_mutation TP53          
## [17] TP53biallelic  trisomy12      vhmut          WBC           
## [21] X11q_biallelic X11q_mono_del  X13q_hom      
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##        CP nsplit rel error xerror  xstd
## 1  0.1961      0     1.000   1.07 0.071
## 2  0.0654      1     0.804   0.87 0.070
## 3  0.0490      4     0.608   0.77 0.069
## 4  0.0392      5     0.559   0.71 0.067
## 5  0.0275      7     0.480   0.70 0.067
## 6  0.0196     12     0.343   0.63 0.065
## 7  0.0147     19     0.206   0.67 0.066
## 8  0.0098     21     0.176   0.64 0.066
## 9  0.0074     29     0.098   0.71 0.067
## 10 0.0065     33     0.069   0.74 0.068
## 11 0.0049     36     0.049   0.74 0.068
## 12 0.0033     40     0.029   0.76 0.069
## 13 0.0020     43     0.020   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.89 0.070
## 3 0.0490      4      0.61   0.77 0.069
## 4 0.0392      5      0.56   0.73 0.068
## 5 0.0275      7      0.48   0.74 0.068
## 6 0.0196     12      0.34   0.62 0.065
## 7 0.0131     19      0.21   0.63 0.065
## 8 0.0098     22      0.17   0.63 0.065
## 9 0.0020     26      0.13   0.68 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           plat           SAMHD1        
##  [9] SF3B1_mutation TP53           WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##        CP nsplit rel error xerror  xstd
## 1  0.1961      0      1.00   1.07 0.071
## 2  0.0654      1      0.80   0.92 0.071
## 3  0.0490      4      0.61   0.88 0.070
## 4  0.0392      5      0.56   0.77 0.069
## 5  0.0275      7      0.48   0.73 0.068
## 6  0.0196     12      0.34   0.70 0.067
## 7  0.0131     14      0.30   0.74 0.068
## 8  0.0098     17      0.26   0.74 0.068
## 9  0.0049     20      0.24   0.81 0.069
## 10 0.0020     22      0.23   0.80 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] trisomy12      vhmut          WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.16 0.070
## 2 0.0654      1      0.80   0.89 0.070
## 3 0.0490      4      0.61   0.81 0.069
## 4 0.0392      5      0.56   0.78 0.069
## 5 0.0245      7      0.48   0.72 0.068
## 6 0.0196      9      0.43   0.69 0.067
## 7 0.0131     14      0.33   0.68 0.067
## 8 0.0049     17      0.29   0.70 0.067
## 9 0.0020     19      0.28   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.91 0.070
## 3 0.0490      4      0.61   0.92 0.071
## 4 0.0392      7      0.46   0.79 0.069
## 5 0.0294      8      0.42   0.75 0.068
## 6 0.0261     10      0.36   0.72 0.068
## 7 0.0049     13      0.28   0.74 0.068
## 8 0.0020     15      0.27   0.74 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.85 0.070
## 3 0.0490      4      0.61   0.81 0.069
## 4 0.0294      7      0.46   0.76 0.069
## 5 0.0261      8      0.43   0.76 0.069
## 6 0.0098     11      0.35   0.77 0.069
## 7 0.0020     15      0.31   0.78 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.92 0.071
## 3 0.0490      4      0.61   0.94 0.071
## 4 0.0294      7      0.46   0.85 0.070
## 5 0.0147      8      0.43   0.81 0.069
## 6 0.0098     10      0.40   0.83 0.070
## 7 0.0020     14      0.36   0.78 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho       ANC_Neutro       CNAs             haem            
## [5] independent_13q. SF3B1_mutation   TP53             WBC             
## [9] X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.88 0.070
## 3 0.0490      4      0.61   0.78 0.069
## 4 0.0392      6      0.51   0.75 0.068
## 5 0.0294      7      0.47   0.75 0.068
## 6 0.0098      8      0.44   0.75 0.068
## 7 0.0065     10      0.42   0.71 0.067
## 8 0.0020     13      0.40   0.71 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.09 0.071
## 2 0.065      1      0.80   0.90 0.070
## 3 0.049      4      0.61   0.83 0.070
## 4 0.025      6      0.51   0.71 0.067
## 5 0.002      8      0.46   0.72 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.91 0.070
## 3 0.049      4      0.61   0.90 0.070
## 4 0.020      6      0.51   0.81 0.069
## 5 0.002      8      0.47   0.80 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.86 0.070
## 3 0.049      4      0.61   0.76 0.069
## 4 0.020      6      0.51   0.76 0.069
## 5 0.002      8      0.47   0.77 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   1.00 0.071
## 3 0.049      4      0.61   0.85 0.070
## 4 0.020      6      0.51   0.93 0.071
## 5 0.002      8      0.47   0.96 0.071
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.065      1      0.80   1.05 0.071
## 3 0.049      4      0.61   0.97 0.071
## 4 0.020      6      0.51   0.88 0.070
## 5 0.002      8      0.47   0.90 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.16 0.070
## 2 0.083      1      0.80   0.93 0.071
## 3 0.039      3      0.64   0.76 0.069
## 4 0.020      6      0.52   0.82 0.069
## 5 0.002      7      0.50   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.083      1      0.80   0.82 0.069
## 3 0.039      3      0.64   0.78 0.069
## 4 0.020      6      0.52   0.81 0.069
## 5 0.002      7      0.50   0.79 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     Binet          BIRC3_mutation
##  [5] cd38           clones         CNAs           del4p         
##  [9] del8p          gender         haem           MED12mutation 
## [13] plat           SAMHD1         SF3B1_mutation TP53          
## [17] TP53biallelic  trisomy12      vhmut          WBC           
## [21] X11q_biallelic X11q_mono_del  X13q_hom      
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##        CP nsplit rel error xerror  xstd
## 1  0.1961      0     1.000   1.07 0.071
## 2  0.0654      1     0.804   0.89 0.070
## 3  0.0490      4     0.608   0.74 0.068
## 4  0.0392      5     0.559   0.75 0.068
## 5  0.0275      7     0.480   0.79 0.069
## 6  0.0196     12     0.343   0.75 0.068
## 7  0.0147     19     0.206   0.69 0.067
## 8  0.0098     21     0.176   0.70 0.067
## 9  0.0074     29     0.098   0.80 0.069
## 10 0.0065     33     0.069   0.84 0.070
## 11 0.0049     36     0.049   0.84 0.070
## 12 0.0033     40     0.029   0.84 0.070
## 13 0.0010     43     0.020   0.84 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           MED12mutation  plat          
##  [9] SAMHD1         TP53           TP53biallelic  WBC           
## [13] X11q_biallelic X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   0.94 0.071
## 3 0.0490      4      0.61   0.90 0.070
## 4 0.0392      5      0.56   0.85 0.070
## 5 0.0275      7      0.48   0.80 0.069
## 6 0.0196     12      0.34   0.76 0.069
## 7 0.0131     19      0.21   0.81 0.069
## 8 0.0098     22      0.17   0.81 0.069
## 9 0.0010     26      0.13   0.81 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     BIRC3_mutation cd38          
##  [5] CNAs           haem           plat           SAMHD1        
##  [9] SF3B1_mutation TP53           WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##        CP nsplit rel error xerror  xstd
## 1  0.1961      0      1.00   1.00 0.071
## 2  0.0654      1      0.80   0.88 0.070
## 3  0.0490      4      0.61   0.84 0.070
## 4  0.0392      5      0.56   0.81 0.069
## 5  0.0275      7      0.48   0.77 0.069
## 6  0.0196     12      0.34   0.72 0.068
## 7  0.0131     14      0.30   0.71 0.067
## 8  0.0098     17      0.26   0.71 0.067
## 9  0.0049     20      0.24   0.69 0.067
## 10 0.0010     22      0.23   0.68 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] trisomy12      vhmut          WBC            X11q_biallelic
## [13] X11q_mono_del 
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.99 0.071
## 3 0.0490      4      0.61   0.88 0.070
## 4 0.0392      5      0.56   0.85 0.070
## 5 0.0245      7      0.48   0.75 0.068
## 6 0.0196      9      0.43   0.72 0.068
## 7 0.0131     14      0.33   0.68 0.067
## 8 0.0049     17      0.29   0.68 0.067
## 9 0.0010     19      0.28   0.69 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho     ANC_Neutro     cd38           CNAs          
##  [5] haem           plat           SAMHD1         TP53          
##  [9] WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.84 0.070
## 3 0.0490      4      0.61   0.83 0.070
## 4 0.0392      7      0.46   0.75 0.068
## 5 0.0294      8      0.42   0.73 0.068
## 6 0.0261     10      0.36   0.69 0.067
## 7 0.0049     13      0.28   0.72 0.068
## 8 0.0010     15      0.27   0.79 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.00 0.071
## 2 0.0654      1      0.80   0.93 0.071
## 3 0.0490      4      0.61   0.78 0.069
## 4 0.0294      7      0.46   0.73 0.068
## 5 0.0261      8      0.43   0.70 0.067
## 6 0.0098     11      0.35   0.72 0.068
## 7 0.0010     15      0.31   0.71 0.067
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
##  [1] ALC_Lympho       ANC_Neutro       clones           CNAs            
##  [5] haem             independent_13q. NOTCH1_mutation  SF3B1_mutation  
##  [9] TP53             WBC              X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.07 0.071
## 2 0.0654      1      0.80   0.86 0.070
## 3 0.0490      4      0.61   0.77 0.069
## 4 0.0294      7      0.46   0.87 0.070
## 5 0.0147      8      0.43   0.81 0.069
## 6 0.0098     10      0.40   0.81 0.069
## 7 0.0010     14      0.36   0.80 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho       ANC_Neutro       CNAs             haem            
## [5] independent_13q. SF3B1_mutation   TP53             WBC             
## [9] X11q_biallelic  
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##       CP nsplit rel error xerror  xstd
## 1 0.1961      0      1.00   1.09 0.071
## 2 0.0654      1      0.80   0.97 0.071
## 3 0.0490      4      0.61   0.91 0.070
## 4 0.0392      6      0.51   0.87 0.070
## 5 0.0294      7      0.47   0.77 0.069
## 6 0.0098      8      0.44   0.79 0.069
## 7 0.0065     10      0.42   0.82 0.069
## 8 0.0010     13      0.40   0.82 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.96 0.071
## 3 0.049      4      0.61   0.79 0.069
## 4 0.025      6      0.51   0.76 0.069
## 5 0.001      8      0.46   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           WBC            X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.11 0.071
## 2 0.065      1      0.80   1.01 0.071
## 3 0.049      4      0.61   0.83 0.070
## 4 0.020      6      0.51   0.77 0.069
## 5 0.001      8      0.47   0.75 0.068
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.07 0.071
## 2 0.065      1      0.80   0.88 0.070
## 3 0.049      4      0.61   0.86 0.070
## 4 0.020      6      0.51   0.79 0.069
## 5 0.001      8      0.47   0.81 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.11 0.071
## 2 0.065      1      0.80   0.87 0.070
## 3 0.049      4      0.61   0.84 0.070
## 4 0.020      6      0.51   0.89 0.070
## 5 0.001      8      0.47   0.87 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] TP53           X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.065      1      0.80   0.92 0.071
## 3 0.049      4      0.61   0.81 0.069
## 4 0.020      6      0.51   0.77 0.069
## 5 0.001      8      0.47   0.76 0.069
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.00 0.071
## 2 0.083      1      0.80   0.89 0.070
## 3 0.039      3      0.64   0.81 0.069
## 4 0.020      6      0.52   0.80 0.069
## 5 0.001      7      0.50   0.83 0.070
## 
## Classification tree:
## rpart(formula = treegenclin$MRD ~ ., data = treegenclin, method = "class", 
##     control = rpart.control(minbucket = i, xval = 10, cp = cp))
## 
## Variables actually used in tree construction:
## [1] ALC_Lympho     ANC_Neutro     CNAs           haem          
## [5] vhmut          X11q_biallelic
## 
## Root node error: 102/209 = 0.49
## 
## n= 209 
## 
##      CP nsplit rel error xerror  xstd
## 1 0.196      0      1.00   1.14 0.070
## 2 0.083      1      0.80   0.81 0.069
## 3 0.039      3      0.64   0.73 0.068
## 4 0.020      6      0.52   0.80 0.069
## 5 0.001      7      0.50   0.79 0.069
```

```r
# dev.off()


# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v3/figures_trees_03062014/02_tree_alldata_cp001_buk1.pdf')
# plot best tree with these variables
rp.treegenclinbuk3cp008 <- rpart(treegenclin$MRD ~ ., method = "class", data = treegenclin, 
    control = rpart.control(minbucket = 1, xval = 10, cp = 0.001))
prp(rp.treegenclinbuk3cp008, extra = 8, uniform = TRUE, branch = 1, left = FALSE, 
    varlen = 0, main = "All data without age, treatment, freqvhmut \n cp=0.001, endgroup=1, MissClassErr~0.01% ", 
    )
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-22.png) 

```r
# dev.off()

# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/03_tree_alldata_cp01_buk5.pdf')
rp.treegenclinbuk5cp01 <- rpart(treegenclin$MRD ~ ., method = "class", data = treegenclin, 
    control = rpart.control(minbucket = 5, xval = 10, cp = 0.01))
prp(rp.treegenclinbuk5cp01, extra = 2, uniform = TRUE, branch = 1, left = FALSE, 
    varlen = 0, main = "All data without age, treatment, freqvhmut \n cp=0.01, endgroup=5, MissClassErr~15% ", 
    )
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-23.png) 

```r
# dev.off()


# rp.treegenclinbuk2cp001 <- rpart(treegenclin$MRD~., method='class',
# data=treegenclin, control=rpart.control(minbucket=2, xval=10, cp=0.001))
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/04_tree_alldata_cp001_buk2.pdf')
# prp(rp.treegenclinbuk3cp008, extra=2, uniform=TRUE, branch=1, left=FALSE,
# varlen=0, main='All data without age, treatment, freqvhmut \n cp=0.008,
# endgroup=3, MissClassErr~21% ', ) dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/05_tree_alldata_cp001_buk2.pdf')
# prp(rp.treegenclinbuk2cp001, extra=8, uniform=TRUE, branch=1, left=FALSE,
# varlen=0, main='All data without age, treatment, freqvhmut \n cp=0.001,
# endgroup=2, MissClassErr~12% ', ) dev.off()
```

Checking frequency distributions of blood data:
===============================================

```r

hist(genclinv3$haem, breaks = 24)
abline(v = 11, col = 3, lty = 1)
abline(v = 9.1, col = 4, lty = 1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r

hist(genclinv3$plat, breaks = 24)
abline(v = 148, col = 4, lty = 1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 

```r

hist(genclinv3$WBC, breaks = 20)
abline(v = 289, col = 4, lty = 1)
abline(v = 76, col = 4, lty = 1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-33.png) 

```r

hist(genclinv3$ANC_Neutro, breaks = 15)
abline(v = 52, col = 4, lty = 1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-34.png) 

```r

hist(genclinv3$ALC_Lympho, breaks = 20)
abline(v = 20, col = 4, lty = 1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-35.png) 

Principal component analysis for blood data:
============================================

Objective:
Check if all blood variables point in the same direction, e.g. could be reduced to one variable


```r

require(FactoMineR)
```

```
## Loading required package: FactoMineR
```

```r
# PCA with function PCA

pdf()
# Classify patients
blood_pca <- PCA(genclinv3[50:54], scale.unit = TRUE, ncp = 5, graph = T)
# scale all the features, ncp: number of dimensions kept in the results (by
# default 5)

blood_load <- dimdesc(blood_pca)
# patient_load
blood_importance_c1 <- as.data.frame(blood_load$Dim.1$quanti)





# William 2011 (Clin. Oncol.)  ========================== - 17q del (FISH) -
# 11q del (FISH) - vhmut - gender - tri12 - cd38 - 12q del (FISH) ```{r}
# treewilliam <- genclinv2[c(6, 8, 9, 10, 19, 46, 35)]
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/06_MissClassErr_william.pdf')
# par(mfrow=c(1, 1)) for (cp in c(0.01, 0.009, 0.008, 0.006, 0.004, 0.002,
# 0.001)){ source('lastfunction.R') nobuckets <- NULL MissClassErr <- NULL
# CorrectClass <- NULL rootnoderr <- NULL relErr <- NULL for (i in 1:20) { #
# built the tree rp.treewilliam <- rpart(treewilliam$MRD~., method='class',
# data=treewilliam, control=rpart.control(minbucket=i, xval=10, cp=cp)) #
# get the number of buckets nobuckets[i] <- i # get the root node error
# rootnoderr[i] <- rp.treewilliam$parms$prior[2] # Get the rel error for
# each tree rpsum <-printcp(rp.treewilliam) relErr[i] <-last(rpsum) #
# Calculate missclassification error MissClassErr[i]
# <-rootnoderr[i]*relErr[i] # Calculate correctly classified patients
# CorrectClass[i] <- 209 - MissClassErr[i]*209 } if (cp==0.01){
# plot(nobuckets, MissClassErr, pch=0, type='o', ylim=c(0, 0.4),
# main='Variables used in William 2011', xlab='No. of patients in end
# group\n(more patients=less complex model)', ylab='Missclassification
# error=fraction of patients that will be put in wrong group ') legend(12,
# 0.2, c('0.01', '0.009', '0.008', '0.006', '0.004', '0.002', '0.001'),
# cex=0.8, col=c('black','green', 'blue', 'orange', 'purple',
# 'darkslategray4', 'red'), pch=0:7, title='complexity:\n lower = more
# complex') } if (cp==0.009){ lines(x=nobuckets+0.05, MissClassErr,
# type='o', pch=1, col='green') } if (cp==0.008){ lines(x=nobuckets+0.1,
# MissClassErr, type='o', pch=2, col='blue') } if (cp==0.006){
# lines(x=nobuckets+0.15, MissClassErr, type='o', pch=3, col='orange') } if
# (cp==0.004){ lines(x=nobuckets+0.20, MissClassErr, type='o', pch=4,
# col='purple') } if (cp==0.002){ lines(x=nobuckets+0.25, MissClassErr,
# type='o', pch=5, col='darkslategray4') } if (cp==0.001){
# lines(x=nobuckets+0.30,MissClassErr, type='o', pch=6, col='red') } }
# dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/07_tree_william_cp001_buk4.pdf')
# #plot best tree with these variables rp.treewilliam <-
# rpart(treewilliam$MRD~., method='class', data=treewilliam,
# control=rpart.control(minbucket=4, xval=10, cp=0.01)) prp(rp.treewilliam,
# extra=8, uniform=TRUE, branch=1, left=FALSE, varlen=0, main='Variables
# used in William 2011 \n cp=0.001, endgroup=4, MissClassErr~34% ')
# dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/08_tree_william_cp001_buk4.pdf')
# prp(rp.treewilliam, extra=2, uniform=TRUE, branch=1, left=FALSE, varlen=0,
# main='Variables used in William 2011 \n cp=0.001, endgroup=4,
# MissClassErr~34% ') dev.off() ``` Rossi 2014 (Blood)
# ========================== - TP53 disruption *(TP53_mutationONLY,
# TP53biallelic, TP53)* - BIRC3 disruption (45) - SF3B1 mutation - NOTCH1
# mutation - 11q del (46) ```{r} treerossi <- genclinv2[c(9, 15, 16, 17, 18,
# 46, 49, 45)]
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/09_MissClassErr_rossi.pdf')
# par(mfrow=c(1, 1)) for (cp in c(0.01, 0.008, 0.006, 0.004, 0.002, 0.001)){
# missclassbybucket <- NULL nobuckets <- NULL MissClassErr <- NULL
# CorrectClass <- NULL rootnoderr <- NULL for (i in 1:20) { # built the tree
# rp.treerossi <- rpart(treerossi$MRD~., method='class', data=treerossi,
# control=rpart.control(minbucket=i, xval=10, cp=cp)) # get the number of
# buckets nobuckets[i] <- i # get the root node error rootnoderr[i] <-
# rp.treerossi$parms$prior[2] # Get the rel error for each tree rpsum
# <-printcp(rp.treerossi) relErr <-last(rpsum) # Calculate
# missclassification error MissClassErr[i] <-rootnoderr[i]*relErr #
# Calculate correctly classified patients CorrectClass[i] <- 209 -
# MissClassErr[i]*209 } if (cp==0.01){ plot(nobuckets, MissClassErr, pch=0,
# type='o', ylim=c(0, 0.6), main='Variables used in Rossi 2014', xlab='No.
# of patients in end group\n(more patients=less complex model)',
# ylab='Missclassification error=fraction of patients that will be put in
# wrong group ') legend(12, 0.3, c('0.01', '0.008', '0.006', '0.004',
# '0.002', '0.001'), cex=0.8,col=c('black','green', 'blue', 'orange',
# 'purple', 'red'), pch=0:5, title='complexity \n low-> more complex') } if
# (cp==0.008){ lines(MissClassErr, type='o', pch=1, col='green') } if
# (cp==0.006){ lines(MissClassErr, type='o', pch=2, col='blue') } if
# (cp==0.004){ lines(MissClassErr, type='o', pch=3, col='orange') } if
# (cp==0.002){ lines(MissClassErr, type='o', pch=4, col='purple') } if
# (cp==0.001){ lines(MissClassErr, type='o', pch=5, col='red') } } dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/10_tree_rossi_cp001_buk2.pdf')
# #plot best tree with these variables rp.treerossi <-
# rpart(treerossi$MRD~., method='class', data=treerossi,
# control=rpart.control(minbucket=2, xval=10, cp=0.001)) prp(rp.treerossi,
# extra=2, uniform=TRUE, branch=1, left=FALSE, varlen=0, main='Variables
# used in Rossi 2014 \n cp=0.001, minbucket= 2, MissClassErr~38%')
# dev.off() ``` Rosenquist 2014 (Leukemia) ========================== - TP53
# *(TP53_mutationONLY, TP53biallelic, TP53)* - BIRC3 - SF3B1 - NOTCH1 - del
# 17p (35)*(X17q11q12)* - del 11q (46) - trisomy 12 (19) - del 13q (39:41)
# *(X13q_hom. X13q_het, X13q_redefinedMRD_loss_4)* - MYD88 (24) - CNA (12)
# ```{r} treerosenq <- genclinv2[c(9, 15, 16, 17, 18, 19, 45, 24, 35, 39:41,
# 46, 49, 12)]
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/11_MissClassErr_rosenquist.pdf')
# par(mfrow=c(1, 1)) for (cp in c(0.01, 0.008, 0.006, 0.004, 0.002, 0.001)){
# missclassbybucket <- NULL nobuckets <- NULL MissClassErr <- NULL
# CorrectClass <- NULL rootnoderr <- NULL for (i in 1:20) { # built the tree
# rp.treerosenq <- rpart(treerosenq$MRD~., method='class', data=treerosenq,
# control=rpart.control(minbucket=i, xval=10, cp=cp)) # get the number of
# buckets nobuckets[i] <- i # get the root node error rootnoderr[i] <-
# rp.treerosenq$parms$prior[2] # Get the rel error for each tree rpsum
# <-printcp(rp.treerosenq) relErr <-last(rpsum) # Calculate
# missclassification error MissClassErr[i] <-rootnoderr[i]*relErr #
# Calculate correctly classified patients CorrectClass[i] <- 209 -
# MissClassErr[i]*209 } if (cp==0.01){ plot(nobuckets, MissClassErr, pch=0,
# type='o', ylim=c(0, 0.5), main='Variables used in Rosenquist 2014',
# xlab='No. of patients in end group\n(more patients=less complex model)',
# ylab='Missclassification error=fraction of patients that will be put in
# wrong group ') legend(12, 0.3, c('0.01', '0.008', '0.006', '0.004',
# '0.002', '0.001'), cex=0.8,col=c('black','green', 'blue', 'orange',
# 'purple', 'red'), pch=0:5, title='complexity \n low-> more complex') } if
# (cp==0.008){ lines(MissClassErr, type='o', pch=1, col='green') } if
# (cp==0.006){ lines(MissClassErr, type='o', pch=2, col='blue') } if
# (cp==0.004){ lines(MissClassErr, type='o', pch=3, col='orange') } if
# (cp==0.002){ lines(MissClassErr, type='o', pch=4, col='purple') } if
# (cp==0.001){ lines(MissClassErr, type='o', pch=5, col='red') } } dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/12_tree_rosenq_cp001_buk1.pdf')
# #plot best tree with these variables rp.treerosenq <-
# rpart(treerosenq$MRD~., method='class', data=treerosenq,
# control=rpart.control(minbucket=1, xval=10, cp=0.001)) prp(rp.treerosenq,
# extra=2, uniform=TRUE, branch=1, left=FALSE, varlen=0, main='Variables
# used in Rosenquist 2014 \n cp=0.001, minbucket= 1, MissClassErr~25%')
# dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/13_tree_rosenq_cp001_buk3.pdf')
# rp.treerosenq <- rpart(treerosenq$MRD~., method='class', data=treerosenq,
# control=rpart.control(minbucket=3, xval=10, cp=0.001)) prp(rp.treerosenq,
# extra=2, uniform=TRUE, branch=1, left=FALSE, varlen=0, main='Variables
# used in Rosenquist 2014 \n cp=0.001, minbucket= 3, MissClassErr~25%')
# dev.off() ``` Health Economic Plan (HEP) ========================== - TP53
# *(TP53_mutationONLY, TP53biallelic, TP53)* - Isolated del_13q (39:41)
# *(X13q_hom. X13q_het, X13q_redefinedMRD_loss_4)* - Isolated trisomy 12
# (19) - NOTCH1 (18) - del_11q (46) - SF3B1 (17) - BIRC3 del
# *(BIRC3_mutation)* - subclones (11) - SAMHD1 (13, 14) - MED12 - TP53 -
# BIRC3 - SF3B1 (17) - NOTCH1 (18) - del 17p (35)*(X17q11q12)* - del 11q
# (46) - trisomy 12 (19) - del 13q - MYD88 (24) ```{r} treehep <-
# genclinv2[c(9, 15, 16, 49, 39:41, 19, 18, 46, 17, 11, 13, 14, 25)]
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/14_MissClassErr_hep.pdf')
# par(mfrow=c(1, 1)) for (cp in c(0.01, 0.008, 0.006, 0.004, 0.002, 0.001)){
# missclassbybucket <- NULL nobuckets <- NULL MissClassErr <- NULL
# CorrectClass <- NULL rootnoderr <- NULL for (i in 1:20) { # built the tree
# rp.treehep <- rpart(treehep$MRD~., method='class', data=treehep,
# control=rpart.control(minbucket=i, xval=10, cp=cp)) # get the number of
# buckets nobuckets[i] <- i # get the root node error rootnoderr[i] <-
# rp.treehep$parms$prior[2] # Get the rel error for each tree rpsum
# <-printcp(rp.treehep) relErr <-last(rpsum) # Calculate missclassification
# error MissClassErr[i] <-rootnoderr[i]*relErr # Calculate correctly
# classified patients CorrectClass[i] <- 209 - MissClassErr[i]*209 } if
# (cp==0.01){ plot(nobuckets, MissClassErr, pch=0, type='o', ylim=c(0, 0.5),
# main='Variables proposed for Health Economic Plan', xlab='No. of patients
# in end group\n(more patients=less complex model)',
# ylab='Missclassification error=fraction of patients that will be put in
# wrong group ') legend(12, 0.3, c('0.01', '0.008', '0.006', '0.004',
# '0.002', '0.001'), cex=0.8,col=c('black','green', 'blue', 'orange',
# 'purple', 'red'), pch=0:5, title='complexity \n low-> more complex') } if
# (cp==0.008){ lines(MissClassErr, type='o', pch=1, col='green') } if
# (cp==0.006){ lines(MissClassErr, type='o', pch=2, col='blue') } if
# (cp==0.004){ lines(MissClassErr, type='o', pch=3, col='orange') } if
# (cp==0.002){ lines(MissClassErr, type='o', pch=4, col='purple') } if
# (cp==0.001){ lines(MissClassErr, type='o', pch=5, col='red') } } dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/15_tree_hep_cp001_buk1.pdf')
# #plot best tree with these variables rp.treehep <- rpart(treehep$MRD~.,
# method='class', data=treehep, control=rpart.control(minbucket=1, xval=10,
# cp=0.001)) prp(rp.treehep, extra=8, uniform=TRUE, branch=1, left=FALSE,
# varlen=0, main='Variables proposed for Health Economics Plan \n cp=0.001,
# endgroup=1, MissClassErr~27%') dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/16_tree_hep_cp001_buk1.pdf')
# prp(rp.treehep, extra=2, uniform=TRUE, branch=1, left=FALSE, varlen=0,
# main='Variables proposed for Health Economics Plan \n cp=0.001,
# endgroup=1, MissClassErr~27%') dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/17_tree_hep_cp001_buk4.pdf')
# rp.treehepb <- rpart(treehep$MRD~., method='class', data=treehep,
# control=rpart.control(minbucket=4, xval=10, cp=0.001)) prp(rp.treehepb,
# extra=2, uniform=TRUE, branch=1, left=FALSE, varlen=0, main='Variables
# proposed for Health Economics Plan \n cp=0.001, endgroup=4,
# MissClassErr~31%') dev.off() ``` Only genetic data =================
# ```{r} treegenetic <- genclinv2[c(9, 11:49)]
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/18_MissClassErr_genetics.pdf')
# par(mfrow=c(1, 1)) for (cp in c(0.01, 0.009, 0.008, 0.007, 0.006, 0.004,
# 0.002, 0.001)){ nobuckets <- NULL MissClassErr <- NULL CorrectClass <-
# NULL rootnoderr <- NULL relErr <- NULL print(cp) for (i in 1:15) { # built
# the tree rp.treegenetic <- rpart(treegenetic$MRD~., method='class',
# data=treegenetic, control=rpart.control(minbucket=i, xval=10, cp=cp)) #
# get the number of buckets nobuckets[i] <- i # get the root node error
# rootnoderr[i] <- rp.treegenetic$parms$prior[2] # Get the rel error for
# each tree rpsum <-printcp(rp.treegenetic) relErr[i] <-last(rpsum) #
# Calculate missclassification error MissClassErr[i]
# <-rootnoderr[i]*relErr[i] # Calculate correctly classified patients
# CorrectClass[i] <- 209 - MissClassErr[i]*209 } if (cp==0.01){
# plot(nobuckets, MissClassErr, pch=0, type='o', ylim=c(0, 0.4), main='Only
# our genetic findings', xlab='No. of patients in end group\n(more
# patients=less complex model)', ylab='Missclassification error=fraction of
# patients that will be put in wrong group') legend(12, 0.25, c('0.01',
# '0.009', '0.008', '0.007', '0.006', '0.004', '0.002', '0.001'), cex=0.8,
# col=c('black', 'green', 'darkslategray4', 'darkgoldenrod', 'blue',
# 'orange', 'purple', 'red'), pch=0:7, title='complexity \n low-> more
# complex') } if (cp==0.009){ lines(MissClassErr, type='o', pch=1,
# col='green') } if (cp==0.008){ lines(MissClassErr, type='o', pch=2,
# col='darkslategray4') } if (cp==0.007){ lines(MissClassErr, type='o',
# pch=3, col='darkgoldenrod') } if (cp==0.006){ lines(MissClassErr,
# type='o', pch=4, col='blue') } if (cp==0.004){ lines(MissClassErr,
# type='o', pch=5, col='orange') } if (cp==0.002){ lines(MissClassErr,
# type='o', pch=6, col='purple') } if (cp==0.001){ lines(MissClassErr,
# type='o', pch=7, col='red') } } dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/19_tree_genetic_cp001_buk1.pdf')
# #plot best tree with these variables rp.treegenetic <-
# rpart(treegenetic$MRD~., method='class', data=treegenetic,
# control=rpart.control(minbucket=1, xval=10, cp=0.001)) prp(rp.treegenetic,
# extra=2, uniform=TRUE, branch=1, left=FALSE, varlen=0, main='Only our
# genetic findings \n cp=0.001, endgroup=1, MissClassErr~11%') dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/20_tree_genetic_cp001_buk3.pdf')
# rp.treegeneticb <- rpart(treegenetic$MRD~., method='class',
# data=treegenetic, control=rpart.control(minbucket=3, xval=10, cp=0.001))
# prp(rp.treegeneticb, extra=2, uniform=TRUE, branch=1, left=FALSE,
# varlen=0, main='Only our genetic findings \n cp=0.001, endgroup=5,
# MissClassErr~22%') dev.off()
# pdf('~/work/01_HICF1/HICF1_sub1/trunk/HICF_predictivemodel/HICF1_trees/HICF1_trees_v2/figures_trees_29052014/21_tree_genetic_cp001_buk8.pdf')
# rp.treegeneticc <- rpart(treegenetic$MRD~., method='class',
# data=treegenetic, control=rpart.control(minbucket=3, xval=10, cp=0.001))
# prp(rp.treegeneticc, extra=2, uniform=TRUE, branch=1, left=FALSE,
# varlen=0, main='Only our genetic findings \n cp=0.001, endgroup=8,
# MissClassErr~25%') dev.off()
```
