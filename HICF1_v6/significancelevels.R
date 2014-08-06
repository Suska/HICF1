#Dr. Susanne Weller
#24/07/2014
#
#This converts p-values into * signs for later use in heatmaps

significancelevels <- function(pvalue){
if(is.na(pvalue)){
    sign <- ""
  }
  else{
  if(pvalue<=0.001){
    sign <- "***"
  }
  if(pvalue<=0.01 & pvalue>0.001){
    sign <- "**"
  }
  if(pvalue<=0.05 & pvalue > 0.01){
    sign <- "*"
  }
  if(pvalue<=0.1 & pvalue > 0.05){
    sign <- "t"
  }
  if(pvalue > 0.01){
    sign <- ""
  }  
return(sign)

  } 
}
