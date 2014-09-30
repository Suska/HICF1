#Dr. Susanne Weller
#23/07/2014
#Calculate odds ratio data frame from data frame

oddsratio <- function(a,b){
  
  if(a>=b){
    return(NA)
  }  
  else{
    if (class(genclinv7.ass[[a]]) == "factor" & class(genclinv7.ass[[b]])=="factor"){
      sum_a <- summary(genclinv7.ass[[a]])
      observedfreq_a <- sum_a[2]/length(which(!is.na(genclinv7.ass[[a]])))
      sum_b <-  summary(genclinv7.ass[[b]])
      observedfreq_b<- sum_b[2]/length(which(!is.na(genclinv7.ass[[b]])))
      #This calculates the expected counts:
      expectedfreq_1 <- observedfreq_a*observedfreq_b
      expectedfreq_0 <- 1-expectedfreq_1
      expectedcount_1 <- expectedfreq_1*250
      expectedcount_0 <- expectedfreq_0*250
      #This caluclates the observed counts:
      table_a_b <- table(genclinv7.ass[[a]],genclinv7.ass[[b]])
      observedcount_1 <- table_a_b[4]
      observedcount_0 <- sum(table_a_b) - table_a_b[4]
      
      oddsratio <- (expectedcount_0*observedcount_1)/(expectedcount_1*observedcount_0)
      

      return(oddsratio)
    }
    else{        
        return(NA)    
    }
  }    
} 