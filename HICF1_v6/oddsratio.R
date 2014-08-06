#Dr. Susanne Weller
#23/07/2014
#Calculate odds ratio data frame from data frame

oddsratio <- function(a,b){
  
  if(a>=b){
    return(NA)
  }  
  else{
    if (class(genclinv6[[a]]) == "factor" & class(genclinv6[[b]])=="factor"){
      sum_a <- summary(genclinv6[[a]])
      observedfreq_a <- sum_a[2]/length(which(!is.na(genclinv6[[a]])))
      sum_b <-  summary(genclinv6[[b]])
      observedfreq_b<- sum_b[2]/length(which(!is.na(genclinv6[[b]])))
      #This calculates the expected counts:
      expectedfreq_1 <- observedfreq_a*observedfreq_b
      expectedfreq_0 <- 1-expectedfreq_1
      expectedcount_1 <- expectedfreq_1*239
      expectedcount_0 <- expectedfreq_0*239
      #This caluclates the observed counts:
      table_a_b <- table(genclinv6[[a]],genclinv6[[b]])
      observedcount_1 <- table_a_b[4]
      observedcount_0 <- sum(table_a_b) - table_a_b[4]
      
      oddsratio <- (expectedcount_0*observedcount_1)/(expectedcount_1*observedcount_0)
      
      # scale oddsratio for plot
      
#       if(oddsratio>25){ 
#          oddsratio=30
#        }
#       if(oddsratio>30 & oddsratio < 35){
#         oddsratio=35
#       }
#       if(oddsratio>35 & oddsratio < 40){
#         oddsratio=40
#       }
#       if(oddsratio > 40){
#         oddsratio=45
#       }
      return(oddsratio)
    }
    else{        
        return(NA)    
    }
  }    
} 