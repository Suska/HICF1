#Dr. Susanne Weller
#28/05/2014

# Calculates missclassification error from recursive partitioning trees for a SET cp and minbucket!)
MissClassErr <- function(data, i, cp){
library(rpart)
rootnode <- function(data, i, cp){
#this builts the rpart model
rp.tree <- rpart(data$MRD~., method="class", data=data, control=rpart.control(minbucket=i, xval=10, cp=cp))
# get the root node error
rootnoderr <- rp.tree$parms$prior[2]
return(rootnoderr)
}

relErr <- function(data, i, cp){
  #this builts the rpart model
  rp.tree <- rpart(data$MRD~., method="class", data=data, control=rpart.control(minbucket=i, xval=cp, cp=cp))
  rpsum <-printcp(rp.tree)
  source("lastfunction.R")
  relErr <-last(rpsum)
  return(relErr)
}
rootnode(data, i, cp)*relErr(data, i, cp)
}