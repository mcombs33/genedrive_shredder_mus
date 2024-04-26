### Title: Script to generate variable parameters for Y-Shredder limited sensitivity analysis using latin hypercube sampling
### Author: Matthew Combs
### Date: 26APR2024
### Notes: output values must be combined with static parameters in correct column order prior to running global sensitivity analysis
########################
library(lhs)

#Limited Y-Shredder sensitivity analysis
A<-randomLHS(1000, 6) #6 variables, 1000 samples
A2<-A
A2[,1]<-floor(qunif(A[,1], min=1, max=5))  #max polyandry (1-4 male mates per female)
A2[,2]<-floor(qunif(A[,2], min=1, max=9)) #max polygyny (1-8 female mates per male)
A2[,3]<-floor(qunif(A[,3], min=1, max=4))  #polyandry rate (low, medium, high)
A2[,4]<-qunif(A[,4], min=0.6, max=3.0)     #dispersal kernel rate (0.6 - 3.0)
A2[,5]<-qunif(A[,5], min=0.005, max=0.25)   #juvenile dispersal percentage
A2[,6]<-qunif(A[,6], min=0.005, max=0.25)   #adult dispersal percentage

colnames(A2)<-c("maxAndry", "maxGyny", "andryRate", "kernRate", "jDispFq", "aDispFq")
write.csv(A2, "YShred_SensAnLimited_LHSparm.csv", row.names = F)
