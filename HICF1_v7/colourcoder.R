#Dr. Susanne Weller
#08/09/2014
#############
#COLOURCODER#
#############

#This function was written to turn odds ratios into a colour scheme to be used with heatmap.2 (gplots package)

# colourcoder <- function(x){
#   if(is.na(x)){
#     return(NA)
#   }
#   if(x<=0.25){
#     return(1)
#   }
#   if((x<=0.5)&(x>0.25)){
#     return(2)
#   }
#   if((x<=0.75)&(x>0.5)){
#     return(3)
#   }
#   if((x<=1)&(x>0.75)){
#     return(4)
#   }
#   if((x<=5)&(x>1)){
#     return(5)
#   }
#   if((x<=10)&(x>5)){
#     return(6)
#   }
#   if((x<=15)&(x>10)){
#     return(7)
#   }
#   if((x<=20)&(x>15)){
#     return(8)
#   }
#   if((x<=30)&(x>20)){
#     return(9)
#   }
# }
library(magrittr)
colourcoder <- function(x){
  ifelse(is.na(x), NA, ifelse(x<=0.25, 1, ifelse((x<=0.5)&(x>0.25),2,ifelse((x<=0.75)&(x>0.5),3,
  ifelse((x<=1)&(x>0.75), 4, ifelse((x<=5)&(x>1), 5, ifelse((x<=10)&(x>5), 6, ifelse((x<=15)&(x>10), 7,
  ifelse((x<=20)&(x>15), 8, ifelse((x<=30)&(x>20), 9, 9))))))))))
}



#   if(){
#     return(4)
#   }
#   if(){
#     return(5)
#   }
#   if((x<=10)&(x>5)){
#     return(6)
#   }
#   if(){
#     return(7)
#   }
#   if(){
#     return(8)
#   }
#   if(){
#     return(9)
#   }
# }