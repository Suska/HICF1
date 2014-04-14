HICF1 - connecting clinical data with WGA data
========================================================
  
  Import clinical data


```r
clinicaldata = read.table("HICF1_clinicaldata.txt", header = TRUE, sep = "\t")
```

```
## Warning: cannot open file 'HICF1_clinicaldata.txt': No such file or
## directory
```

```
## Error: cannot open the connection
```

```r
clinicaldata[clinicaldata == "9876"] <- NA
```

```
## Error: object 'clinicaldata' not found
```


Convert date info in format 'mm/dd/yyyy' and remove corresponding string afterwards, recalculate age at randomisation to one digit, remove double columns, rearrange columns


```r
strDates <- clinicaldata$patient_date_of_birth
```

```
## Error: object 'clinicaldata' not found
```

```r
clinicaldata$dob <- as.Date(strDates, "%d/%m/%Y")
```

```
## Error: object 'strDates' not found
```

```r
clinicaldata$patient_date_of_birth <- NULL
```

```
## Error: object 'clinicaldata' not found
```

```r

strDates_dor <- clinicaldata$date_of_randomisation
```

```
## Error: object 'clinicaldata' not found
```

```r
clinicaldata$dor <- as.Date(strDates_dor, "%d/%m/%Y")
```

```
## Error: object 'strDates_dor' not found
```

```r
clinicaldata$date_of_randomisation <- NULL
```

```
## Error: object 'clinicaldata' not found
```

```r

clinicaldata$age_ar <- round(((difftime(clinicaldata$dor, clinicaldata$dob, 
    units = "weeks"))/52), digits = 1)
```

```
## Error: object 'clinicaldata' not found
```

```r
clinicaldata$age_at_randomisation <- NULL
```

```
## Error: object 'clinicaldata' not found
```

```r

clinicaldata$ID <- paste(clinicaldata$trial, clinicaldata$clinical_trial_id)
```

```
## Error: object 'clinicaldata' not found
```

```r
clinicaldata$trial <- NULL
```

```
## Error: object 'clinicaldata' not found
```

```r
clinicaldata$clinical_trial_id <- NULL
```

```
## Error: object 'clinicaldata' not found
```

```r

clinicaldata$gender <- factor(clinicaldata$gender_id, labels = c("male", "female"))
```

```
## Error: object 'clinicaldata' not found
```

```r
clinicaldata$gender_id <- NULL
```

```
## Error: object 'clinicaldata' not found
```

```r

a
```

```
## Error: object 'a' not found
```


Rearrange collumns

Overview of the data
-> possible confounding factors (age, sex, ethnicity)
-> treatment groups
