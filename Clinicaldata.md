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









