### Title: Script to generate variable parameters for X-Shredder limited sensitivity analysis using latin hypercube sampling
### Author: Matthew Combs
### Date: 26APR2024
### Notes: output values must be combined with static parameters in correct column order prior to running limited sensitivity analysis
########################
library(lhs)

#Limited X-Shredder sensitivity analysis
A<-randomLHS(1000, 7) #7 variables, 1000 samples
A2<-A

A2[,1]<-floor(qunif(A[,1], min=1, max=5))  #max polyandry (1-4 male mates per female)
A2[,2]<-floor(qunif(A[,2], min=1, max=9))  #max polygyny (1-8 female mates per male)
A2[,3]<-floor(qunif(A[,3], min=1, max=4))  #polyandry rate (low, medium, high)
A2[,4]<-qunif(A[,4], min=0.6, max=3.0)     #dispersal kernel rate (0.6 - 3.0)
A2[,5]<-qunif(A[,5], min=0.05, max=0.25)   #juvenile dispersal percentage
A2[,6]<-qunif(A[,6], min=0.05, max=0.25)   #adult dispersal percentage
A2[,7]<-qunif(A[,7], min=0.5, max=1.0)     #sperm competitiveness of gd males

colnames(A2)<-c("maxAndry", "maxGyny", "andryRate", "kernRate", "jDispFq","aDispFq", "spermComp")
write.csv(A2, "XShred_SensAnLimited_LHSparm.csv", row.names = F)
