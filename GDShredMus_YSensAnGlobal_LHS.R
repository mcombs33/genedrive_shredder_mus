### Title: Script to generate variable parameters for X-Shredder global sensitivity analysis using latin hypercube sampling
### Author: Matthew Combs
### Date: December2025
### Notes: output values must be combined with static parameters in correct column order prior to running global sensitivity analysis
########################
library(lhs)

#Expanded sensitivity analysis
B<-randomLHS(3000, 15) #15 variables, 2000 samples
B[,1]<-qunif(B[,1], min=1, max=1)       #female mating rate (6-10 litters/year)
B[,2]<-floor(qunif(B[,2], min=10, max=101))   #patch carrying capcity (starting density)
B[,3]<-floor(qunif(B[,3], min=1, max=5))      #max polyandry (1-4 male mates per female)
B[,4]<-floor(qunif(B[,4], min=1, max=9))      #max polygyny (1-10 female mates per male)
B[,5]<-floor(qunif(B[,5], min=1, max=4))      #polyandry rate (low, medium, high)
B[,6]<-qunif(B[,6], min=0.6, max=3.0)         #dispersal kernel rate (0.6 - 2.2)
B[,7]<-qunif(B[,7], min=0.001, max=1)       #juvenile dispersal percentage
B[,8]<-qunif(B[,8], min=0.001, max=1)       #adult dispersal percentage
B[,9]<-floor(qunif(B[,9], min=10, max=1001))  #GD release size (10 - 5000)
B[,10]<-floor(qunif(B[,10], min=1, max=12))   #multiple supp schedule
B[,11]<-floor(qunif(B[,11], min=1, max=3))    #central or random placement
B[,12]<-qunif(B[,12], min=0.9, max=0.9)       #prob of cutting homologous site
B[,13]<-qunif(B[,13], min=0.0, max=0.1)       #prob of NHEJ
B[,14]<-qunif(B[,14], min=0.9, max=0.9)       #prob of shredding Y-chromosome
B[,15]<-qunif(B[,15], min=0.5, max=1.0)       #prob of XO females passing X


colnames(B)<-c("femMate","carCap","maxAndry", "maxGyny", "andryRate", "kernRate", 
                "jDispFq","aDispFq","gdIntroSize", "gdSupp", "gdPlace","cutProb", "nhejProb", "shredProb", "Xpass")
write.csv(B, "YShred_SensAnGlobal_LHSparm", row.names = F)
