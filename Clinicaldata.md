HICF1 - connecting clinical data with WGA data
========================================================

Import clinical data


```r
clinicaldata = read.table("/home/suska/work/HICF1/HICF1_clinicaldata/HICF1_clinicaldata_raw/HICF1_clinicaldata.txt", 
    header = TRUE, sep = "\t")
clinicaldata[clinicaldata == "9876"] <- NA
```


Data cleaning
-------------
- Convert date info in format 'mm/dd/yyyy' & remove corresponding string afterwards
- Recalculate age at randomisation to one digit
- Convert int to factors with the correct levels
- Rearrange columns


```r
strDates <- clinicaldata$patient_date_of_birth
clinicaldata$dob <- as.Date(strDates, "%d/%m/%Y")
clinicaldata$patient_date_of_birth <- NULL

strDates_dor <- clinicaldata$date_of_randomisation
clinicaldata$dor <- as.Date(strDates_dor, "%d/%m/%Y")
clinicaldata$date_of_randomisation <- NULL

clinicaldata$age_ar <- as.numeric(round(((difftime(clinicaldata$dor, clinicaldata$dob, 
    units = "weeks"))/52), digits = 1))
clinicaldata$age_at_randomisation <- NULL

clinicaldata$ID <- paste(clinicaldata$trial, clinicaldata$clinical_trial_id)
clinicaldata$trial <- NULL
clinicaldata$clinical_trial_id <- NULL

clinicaldata$gender <- factor(clinicaldata$gender_id, labels = c("male", "female"))
clinicaldata$gender_id <- NULL

clinicaldata$treatment <- factor(clinicaldata$allocation_id, labels = c("FCR", 
    "FCMR", "FCminiR/FCR"))
clinicaldata$allocation_id <- NULL

clinicaldata$vhmut <- factor(clinicaldata$vh_id, labels = c("mutated", "unmutated", 
    "biclonal", "equivocal"))
clinicaldata$vh_id <- NULL

clinicaldata$vh3 <- factor(clinicaldata$vh3.21, labels = c("mutated", "unmutated"))

clinicaldata$vh_freq <- clinicaldata$vh_mutation_.
clinicaldata$vh_mutation_. <- NULL

clinicaldata$X13q14del <- factor(clinicaldata$X13q14_deleted, labels = c("X13_deleted", 
    "X13_not_deleted"))
clinicaldata$X13q14_deleted <- NULL

clinicaldata$X13q14_freq <- clinicaldata$X13q14_._deleted
clinicaldata$X13q14_._deleted <- NULL

clinicaldata$tri12 <- factor(clinicaldata$trisomy12, labels = c("tri12_present", 
    "tri12_absent"))
clinicaldata$trisomy12 <- NULL

clinicaldata$tri12_freq <- clinicaldata$trisomy_._abnormal
clinicaldata$trisomy_._abnormal <- NULL

clinicaldata$X11q23del <- factor(clinicaldata$X11q23_deleted, labels = c("X11_deleted", 
    "X11_not_deleted"))
clinicaldata$X11q23_deleted <- NULL

clinicaldata$X11q23_freq <- clinicaldata$X11q23_._deleted
clinicaldata$X11q23_._deleted <- NULL

clinicaldata$X17pdel <- factor(clinicaldata$X17p_deleted, labels = c("X17_deleted", 
    "X17_not_deleted"))
clinicaldata$X17p_deleted <- NULL

clinicaldata$X17p_freq <- clinicaldata$X17p_._deleted
clinicaldata$X17p_._deleted <- NULL

clinicaldata$cd38 <- factor(clinicaldata$cd38..., labels = c("positive", "negative"))
clinicaldata$cd38... <- NULL

clinicaldata$mrd_id[clinicaldata$mrd_id == 3] <- NA
clinicaldata$mrd_id[clinicaldata$mrd_id == 4] <- NA
clinicaldata$mrd_id[clinicaldata$mrd_id == 5] <- NA
clinicaldata$mrd_id[clinicaldata$mrd_id == 6] <- NA
clinicaldata$mrd <- factor(clinicaldata$mrd_id, labels = c("positive", "negative"))
clinicaldata$mrd_id <- NULL

clinicaldata$final_response_id[clinicaldata$final_response_id == 6] <- NA
clinicaldata$final_response_id[clinicaldata$final_response_id == 7] <- NA
clinicaldata$final_response_id[clinicaldata$final_response_id == 8] <- NA
clinicaldata$final_response_id[clinicaldata$final_response_id == 9] <- NA
clinicaldata$response <- factor(clinicaldata$final_response_id, labels = c("CR", 
    "CR_incomplete", "PR", "SD", "PD", "CR_Cri"))
clinicaldata$final_response_id <- NULL

clinicaldata <- clinicaldata[c("ID", "gender", "dob", "dor", "age_ar", "treatment", 
    "mrd", "response", "vhmut", "vh3", "vh3.21", "vh_freq", "X13q14del", "X13q14_freq", 
    "tri12", "tri12_freq", "X11q23del", "X11q23_freq", "X17pdel", "X17p_freq", 
    "cd38")]

library(gplots)
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
boxplot2(age_ar ~ gender, data = clinicaldata, main = "Gender ~ age", col = c("blue", 
    "red"), xlab = "gender", ylab = "age", top = TRUE, varwidth = "T")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 
