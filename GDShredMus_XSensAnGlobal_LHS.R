### Title: Script to generate variable parameters for X-Shredder global sensitivity analysis using latin hypercube sampling
### Author: Matthew Combs
### Date: December2025
### Notes: output values must be combined with static parameters in correct column order prior to running global sensitivity analysis
########################

library(lhs)
B2<-randomLHS(3000, 12) #13 variables, 2000 samples
B2[,1]<-floor(qunif(B2[,1], min=20, max=101))   #patch carrying capcity (starting density)
B2[,2]<-floor(qunif(B2[,2], min=1, max=5))      #max polyandry (1-4 male mates per female)
B2[,3]<-floor(qunif(B2[,3], min=1, max=9))      #max polygyny (1-10 female mates per male)
B2[,4]<-floor(qunif(B2[,4], min=1, max=4))      #polyandry rate (low, medium, high)
B2[,5]<-qunif(B2[,5], min=0.6, max=3.0)         #dispersal kernel rate (0.6 - 2.2)
B2[,6]<-qunif(B2[,6], min=0.001, max=1)         #juvenile dispersal percentage
B2[,7]<-qunif(B2[,7], min=0.001, max=1)         #adult dispersal percentage
B2[,8]<-floor(qunif(B2[,8], min=10, max=1001))  #GD release size (10 - 5000)
B2[,9]<-floor(qunif(B2[,9], min=1, max=12))   #multiple supp schedule
B2[,10]<-floor(qunif(B2[,10], min=1, max=3))    #central or random placement
B2[,11]<-qunif(B2[,11], min=0.95, max=95)       #prob of cutting X-chromosome
B2[,12]<-qunif(B2[,12], min=0.5, max=1.0)       #sperm competitiveness of gd males

colnames(B2)<-c("carCap","maxAndry", "maxGyny", "andryRate", "kernRate", 
                "jDispFq","aDispFq","gdIntroSize", "gdSupp", "gdPlace", "shredEff","spermComp")
write.csv(B2, "XShred_SensAnGlobal_LHSparm", row.names = F)

