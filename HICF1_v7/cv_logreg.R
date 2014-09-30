#Dr. Susanne Weller
#27/08/2014

#This function takes a dataframe  ith binary response as input and
#computes repeated sub-sampling missclassification error.

cv_logreg <- function(dataset){
truenegative.test <- NULL
falsenegative.test <- NULL
truepositive.test <- NULL
falsepositive.test <-  NULL
missclasserror.test <- NULL
model.no <- NULL
for (i in (1:1000)){
 model.no[i] <- i
  set.seed(i)
  #Use the sample function to generate a random list of numbers 1:209
  randomise.rows <- sample(1:nrow(dataset))
  divider <- round(nrow(dataset)/10)
  trainindex <- randomise.rows[1:(nrow(dataset)-divider)]
  testindex <- randomise.rows[1:divider]
  train <-genclinv6_genetic3[row.names(genclinv6_genetic3) %in% trainindex,]
  test <- genclinv6_genetic3[row.names(genclinv6_genetic3) %in% testindex,]
  #train the model with the training data set
  train.fit.logreg <-glm(MRD ~ ., family=binomial(logit), data=train)
  
  #predict classes with this model and glue it to test data frame
  test$predprob <- predict(train.fit.logreg, test, type="response")
  test$predclass <- ifelse(test$predprob>0.5, "MRD positive", "MRD negative")
  test$correctclassed <- test$predclass==test$MRD
  
   testclass <- as.data.frame(table(test$correctclassed, test$MRD))
   
   #this gives the true values in the test set:
   #summary.test <- summary(test$MRD)
   #this assigns them to variables:
   #allnegative <- summary.test[2]#these are all MRD+ in the test set
   #allpositive <-summary.test[1]#these are all MRD- in the test set
   #calculate the Missclassification Error from test$MRD and rf.predict:
   truenegative.test[i] <- testclass$Freq[4]/divider
   falsenegative.test[i] <- testclass$Freq[1]/divider
   truepositive.test[i] <- testclass$Freq[2]/divider
   falsepositive.test[i] <-  testclass$Freq[3]/divider
   missclasserror.test[i] <- (testclass$Freq[3]+testclass$Freq[1])/divider
 }
   modelperformance <- as.data.frame(cbind(model.no, truenegative.test, falsenegative.test,
                                                   truepositive.test, falsepositive.test, missclasserror.test))
  return (modelperformance)
}
