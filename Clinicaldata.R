HICF1 - connecting clinical data with WGA data
========================================================
  
  Import clinical data

```{r}

clinicaldata = read.table("HICF1_clinicaldata.txt", header = TRUE, sep = "\t")
clinicaldata[clinicaldata=="9876"]<-NA
```

Convert date info in format 'mm/dd/yyyy' and remove corresponding string afterwards, recalculate age at randomisation to one digit, remove double columns, rearrange columns

```{r}

strDates <- clinicaldata$patient_date_of_birth
clinicaldata$dob <- as.Date(strDates, "%d/%m/%Y")
clinicaldata$patient_date_of_birth <- NULL


strDates_dor <- clinicaldata$date_of_randomisation
clinicaldata$dor <- as.Date(strDates_dor, "%d/%m/%Y")
clinicaldata$date_of_randomisation <- NULL 

clinicaldata$age_ar <-round(((difftime(clinicaldata$dor, clinicaldata$dob, units='weeks'))/ 52), digits = 1)
clinicaldata$age_at_randomisation <- NULL


```

Rearrange collumns

Overview of the data
-> possible confounding factors (age, sex, ethnicity)
-> treatment groups
