### Title: Script to generate variable parameters for X-Shredder global sensitivity analysis using latin hypercube sampling
### Author: Matthew Combs
### Date: 26APR2024
### Notes: output values must be combined with static parameters in correct column order prior to running global sensitivity analysis
########################
library(lhs)

#Global X-Shredder sensitivity analysis
B<-randomLHS(2500, 13) #13 variables, 2500 samples
B2<-B
B2[,1]<-qunif(B[,1], min=0.33, max=0.55)       #female mating rate (6-10 litters/year)
B2[,2]<-floor(qunif(B[,2], min=10, max=101))   #patch carrying capcity (starting density)
B2[,3]<-floor(qunif(B[,3], min=1, max=5))      #max polyandry (1-4 male mates per female)
B2[,4]<-floor(qunif(B[,4], min=1, max=9))      #max polygyny (1-8 female mates per male)
B2[,5]<-floor(qunif(B[,5], min=1, max=4))      #polyandry rate (low, medium, high)
B2[,6]<-qunif(B[,6], min=0.6, max=3.0)         #dispersal kernel rate (0.6 - 3)
B2[,7]<-qunif(B[,7], min=0.05, max=0.25)       #juvenile dispersal percentage
B2[,8]<-qunif(B[,8], min=0.05, max=0.25)       #adult dispersal percentage
B2[,9]<-floor(qunif(B[,9], min=10, max=1001))  #GD release size (10 - 1000)
B2[,10]<-floor(qunif(B[,10], min=1, max=12))   #multiple supp schedule
B2[,11]<-floor(qunif(B[,11], min=1, max=3))    #central or random placement
B2[,12]<-qunif(B[,12], min=0.9, max=1.0)       #prob of cutting X-chromosome
B2[,13]<-qunif(B[,13], min=0.5, max=1.0)       #sperm competitiveness of gd males

colnames(B2)<-c("femMate","carCap","maxAndry", "maxGyny", "andryRate", "kernRate", 
                "jDispFq","aDispFq","gdIntroSize", "gdSupp", "gdPlace","shredEff", "spermComp")
write.csv(B2, "XShred_SensAnGlobal_LHSparm.csv", row.names = F)