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

clinicaldata$age_ar <- round(((difftime(clinicaldata$dor, clinicaldata$dob, 
    units = "weeks"))/52), digits = 1)
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

clinicaldata$mrd <- factor(clinicaldata$mrd_id, labels = c("positive", "negative", 
    "NA", "NA", "NA"))
```

```
## Warning: duplicated levels in factors are deprecated
```

```r
clinicaldata$mrd_id <- NULL

clinicaldata$response <- factor(clinicaldata$final_response_id, labels = c("CR", 
    "CR_incomplete", "PR", "SD", "PD", "NA", "NA", "NA", "CR_Cri"))
```

```
## Warning: duplicated levels in factors are deprecated
```

```r
clinicaldata$final_response_id <- NULL

```


Rearrange collumns

Overview of the data
-> possible confounding factors (age, sex, ethnicity)
-> treatment groups
