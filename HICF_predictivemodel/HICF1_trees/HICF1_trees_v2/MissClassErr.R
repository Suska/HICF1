#Dr. Susanne Weller
#28/05/2014

# Calculates missclassification error from recursive partitioning trees for a SET cp and minbucket!)
MissClassErr <- function(data, index){
library(rpart)
rootnode <- function(data){
#this builts the rpart model
rp.tree <- rpart(data$MRD~., method="class", data=data, control=rpart.control(minbucket=1, xval=10, cp=0.001))
# get the root node error
rootnoderr <- rp.tree$parms$prior[2]
return(rootnoderr)
}

relErr <- function(data){
  #this builts the rpart model
  rp.tree <- rpart(data$MRD~., method="class", data=data, control=rpart.control(minbucket=1, xval=10, cp=0.001))
  rpsum <-printcp(rp.tree)
  source("lastfunction.R")
  relErr <-last(rpsum)
  return(relErr)
}
rootnode(data[index,])*relErr(data[index,])
}