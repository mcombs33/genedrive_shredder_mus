### Title: Script to run global sensitivity analysis of Y-Shredder gene drive
### Author: Matthew Combs
### Date: 26APR2024
### Notes: parameters established beforehand and uploaded as .csv file
########################
library(parallel)
library(doParallel)
library(deSolve)   
library(data.table)
library(extraDistr)
library(bigstatsr)
library(bigreadr)

#############################################################################
#load functions
source(file = "GDShredMus_Functions.R")

##############################################################################
#Create output directory for csv outputs
outfolder<-"YShred_SensGlobal"
basepath.out<-"/lustrefs/nwrc/projects/gdsim/outfiles/"
outloc<-paste0(basepath.out, outfolder)
dir.create(outloc)

##############################################################################
#Set lattice size, time steps, iterations
NL<-31
weeks<- 742 #30 years + 22 timesteps (burn-in)
iterations<-50

##############################################################################
#load params from csv
tparms.csv<-read.csv(file="YSensAnGlobal_parms.csv", header=T)#import csv of params
tparms.df<-as.data.frame(sapply(tparms.csv, as.numeric)) #make numeric

vec1<-c(1:2500)
tparms.df[1:nrow(tparms.df),(ncol(tparms.df)+1)]<-vec1  

tparms.vec<-as.vector(t(tparms.df)) #transpose and vectorize
tparms.mat<-if(nrow(tparms.csv)==1){
  matrix(data=tparms.vec, nrow=nrow(tparms.df), ncol=ncol(tparms.df), byrow=F) #reform into transposed matrix
} else {
  matrix(data=tparms.vec, nrow=ncol(tparms.df), ncol=nrow(tparms.df), byrow=F) #reform into transposed matrix
}
tparms<-matrix(data=apply(tparms.mat, 2, 
                          function(x) rep(x, iterations)), 
               ncol= ncol(tparms.mat) * iterations) #replicate each parm set by number of iterations

#Setup matrix to write basic results
runResults<-FBM(ncol(tparms),5)

################################################################################
#Establish additional run parameters

#Null mating matrix, used in randBirthYS function
b0 <-rep.int(0, NL*NL*10)
Birthnull<-matrix(b0,ncol=(10),nrow=NL*NL,byrow = T) 

#Mating System setup
AndryProbs1<-c(1) 
AndryProbs2<-c(0.9,0.1) 
AndryProbs3<-c(0.9,0.05,0.05) 
AndryProbs4<-c(0.9,0.05,0.03,0.02) 
AndryProbs5<-c(1) 
AndryProbs6<-c(0.7,0.3) 
AndryProbs7<-c(0.7,0.2,0.1) 
AndryProbs8<-c(0.7,0.15,0.1,0.05) 
AndryProbs9<-c(1) 
AndryProbs10<-c(0.5,0.5) 
AndryProbs11<-c(0.5,0.3,0.2) 
AndryProbs12<-c(0.5,0.25,0.15,0.1) 
PolyandryProbsList<-list(AndryProbs1,AndryProbs2,AndryProbs3,AndryProbs4,
                         AndryProbs5,AndryProbs6,AndryProbs7,AndryProbs8,
                         AndryProbs9,AndryProbs10,AndryProbs11,AndryProbs12)

################################################################################
#Setup cluster
n.cores <-96  #find the number of cores in system
my.cluster <- parallel::makeCluster(n.cores, type = "FORK", outfile="debug.txt")
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

################################################################################
#Run the model
start.time <- Sys.time() #Just to measure how long the function takes.

foreach(
  i = 1:ncol(tparms),
  .combine = 'rbind'
) %dopar% {
  library(deSolve)
  #Establish initial abundances 
  #Adults
  XXw0 <-matrix((floor(tparms[5,i] * 0.5)),NL,NL)   #wt females
  XYw0 <-matrix((floor(tparms[5,i] * 0.5)),NL,NL)   #wt males
  XOw0 <-matrix(0,NL,NL)    #wt XO females
  XXH0 <-matrix(0,NL,NL)    #Hom XX females
  XXh0 <-matrix(0,NL,NL)    #Het XX females
  XYH0 <-matrix(0,NL,NL)    #Hom XY males
  XYh0 <-matrix(0,NL,NL)    #Het XY females
  XOH0 <-matrix(0,NL,NL)    #Hom XO females
  XOH0 <-matrix(0,NL,NL)    #Het XO females
  #Juveniles
  jXXw0 <-matrix(0,NL,NL)   #wt females
  jXYw0 <-matrix(0,NL,NL)   #wt males
  jXOw0 <-matrix(0,NL,NL)   #wt XO females
  jXXH0 <-matrix(0,NL,NL)   #Hom XX females
  jXXh0 <-matrix(0,NL,NL)   #Het XX females
  jXYH0 <-matrix(0,NL,NL)   #Hom XY males
  jXYh0 <-matrix(0,NL,NL)   #Het XY females
  jXOH0 <-matrix(0,NL,NL)   #Hom XO females
  jXOh0 <-matrix(0,NL,NL)   #Het XO females
  
  #Calculate summary values
  Total0 <- XXw0+XYw0+XOw0+XXH0+XXh0+XYH0+XYh0+XOH0+XOH0+jXXw0+jXYw0+jXOw0+jXXH0+jXXh0+jXYH0+jXYh0+jXOH0+jXOh0
  GDMO0 <- XXH0+XXh0+XYH0+XYh0+XOH0+XOH0+jXXH0+jXXh0+jXYH0+jXYh0+jXOH0+jXOh0
  #SexRatio
  sexRatio<-(XYw0+XYH0+XYh0+jXYw0+jXYH0+jXYh0)/Total0
  #Empty Births
  births<-c(rep.int(0,NL*NL))
  #Average percentage of indivduals that are GDMO across patches
  Out10 <- sum(GDMO0)/sum(Total0) 
  #Average percentage of resistant males across patches
  #Out20 <-sum(RY0 + RD0)/sum(Total0)
  Out20<-0
  #Overall abundance of all individuals
  Out30 <- sum(Total0)
  #Combine into list of matrices to initialize model
  yn <-c(XXw0,XYw0,XOw0,XXH0,XXh0,XYH0,XYh0,XOH0,XOH0,jXXw0,jXYw0,jXOw0,jXXH0,jXXh0,jXYH0,jXYh0,jXOH0,jXOh0,
         Total0, sexRatio, births, GDMO0, 
         0,0,Out10,Out20,Out30,0,0,1)
  
  #Establish mating system parameters for specific iteration
  PolyandrySeq<-seq(1,tparms[6,i],1)   #setup for polyandry 
  AndryRate<-tparms[21,i]
  if(AndryRate==1){
    PolyandryProbs<-PolyandryProbsList[[(tparms[6,i])]]
  } else if(AndryRate==2){
    PolyandryProbs<-PolyandryProbsList[[ 4 + (tparms[6,i]) ]]
  } else if(AndryRate==3){
    PolyandryProbs<-PolyandryProbsList[[ 8 + (tparms[6,i]) ]]
  }
  ##############################################################################
  #Setup run specific genotype probabilities for each iteration
  Pc<-tparms[15,i]
  Pn<-tparms[16,i]
  Py<-tparms[17,i]
  Xpass<-tparms[18,i]
  
  #P1: XXw * XYw
  P1XXw<- 0.5
  P1XYw<- 0.5
  P1XXH<- 0
  P1XXh<- 0
  P1XYH<- 0
  P1XYh<- 0
  P1XOH<- 0
  P1XOh<- 0
  P1XOw<- 0
  P1nv <- 0
  
  #P2: XXw * XYh
  P2XXw<- (1-Pc)*0.5*0.5
  P2XYw<- (1-Pc)*0.5*0.5
  P2XXH<- 0
  P2XXh<- (1-Pc)*0.5*0.5 + Pc*(1-Pn)*0.5 + Pc*Pn*0.25
  P2XYH<- 0
  P2XYh<- (1-Pc)*0.5*0.5*(1-Py) + Pc*(1-Pn)*0.5*(1-Py) + Pc*Pn*0.25*(1-Py)
  P2XOH<- 0
  P2XOh<- (1-Pc)*0.5*0.5*Py + Pc*(1-Pn)*0.5*Py +Pc*Pn*0.25*Py
  P2XOw<- 0
  P2nv <- 2*(Pc*Pn*0.25)
  
  #P3: XXw * XYH
  P3XXw<- 0
  P3XYw<- 0
  P3XXH<- 0
  P3XXh<- 0.5
  P3XYH<- 0
  P3XYh<- 0.5*(1-Py)
  P3XOH<- 0
  P3XOh<- 0.5*Py
  P3XOw<- 0
  P3nv <- 0
  
  #P4: XXh * XYw
  P4XXw<- (1-Pc)*0.5*0.5
  P4XYw<- (1-Pc)*0.5*0.5
  P4XXH<- 0
  P4XXh<- (1-Pc)*0.5*0.5 + Pc*(1-Pn)*0.5 + Pc*Pn*0.25
  P4XYH<- 0
  P4XYh<- (1-Pc)*0.5*0.5*(1-Py) + Pc*(1-Pn)*0.5*(1-Py) + Pc*Pn*0.25*(1-Py)
  P4XOH<- 0
  P4XOh<- (1-Pc)*0.5*0.5*Py + Pc*(1-Pn)*0.5*Py + Pc*Pn*0.25*Py
  P4XOw<- 0
  P4nv <- 2*(Pc*Pn*0.25)
  
  #P5: XXh * XYh 
  P5XXw<- (1-Pc)^2*0.125
  P5XYw<- (1-Pc)^2*0.125
  P5XXH<- (1-Pc)*0.5*Pc*(1-Pn)*0.5 + Pc*(1-Pn)*(1-Pc)*0.5*0.5 + (1-Pc)*0.5*(1-Pc)*0.5*0.5 +  Pc*(1-Pn)*Pc*(1-Pn)*0.5 +   Pc*Pn*0.5*Pc*(1-Pn)*0.5 +  Pc*Pn*0.5*(1-Pc)*0.25 +  Pc*(1-Pn)*Pc*Pn*0.25 +  (1-Pc)*0.5*Pc*Pn*0.25 +Pc*Pn*0.5*Pc*Pn*0.25
  P5XXh<- (1-Pc)*0.5*Pc*(1-Pn)*0.5 + Pc*(1-Pn)*(1-Pc)*0.5*0.5 + 2*((1-Pc)*0.5*(1-Pc)*0.5*0.5) +  (1-Pc)*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*(1-Pc)*0.25 
  P5XYH<- (1-Pc)*0.5*Pc*(1-Pn)*0.5*(1-Py)^2 + Pc*(1-Pn)*(1-Pc)*0.5*0.5*(1-Py)^2 + (1-Pc)*0.5*(1-Pc)*0.5*0.5*(1-Py)^2 + Pc*(1-Pn)*Pc*(1-Pn)*0.5*(1-Py)^2 +  Pc*Pn*0.5*Pc*(1-Pn)*0.5*(1-Py)^2 +  Pc*Pn*0.5*(1-Pc)*0.25*(1-Py)^2 +  Pc*(1-Pn)*Pc*Pn*0.25*(1-Py)^2 +  (1-Pc)*0.5*Pc*Pn*0.25*(1-Py)^2 +Pc*Pn*0.5*Pc*Pn*0.25*(1-Py)^2
  P5XYh<- (1-Pc)*0.5*Pc*(1-Pn)*0.5*(1-Py) + Pc*(1-Pn)*(1-Pc)*0.5*0.5*(1-Py) + 2*((1-Pc)*0.5*(1-Pc)*0.5*0.5*(1-Py)) +  (1-Pc)*0.5*Pc*Pn*0.25*(1-Py) +  Pc*Pn*0.5*(1-Pc)*0.25*(1-Py) 
  P5XOH<-   (1-Pc)*0.5*Pc*(1-Pn)*0.5*(1-(1-Py)^2) +   Pc*(1-Pn)*(1-Pc)*0.5*0.5*(1-(1-Py)^2) +   (1-Pc)*0.5*(1-Pc)*0.5*0.5*(1-(1-Py)^2) +   Pc*(1-Pn)*Pc*(1-Pn)*0.5*(1-(1-Py)^2) +  (1-Pc)*0.5*Pc*Pn*0.25*(1-(1-Py)^2) +  Pc*Pn*0.5*Pc*(1-Pn)*0.5*(1-(1-Py)^2) +  Pc*(1-Pn)*Pc*Pn*0.25 *(1-(1-Py)^2) +  Pc*Pn*0.5*(1-Pc)*0.25*(1-(1-Py)^2) +  Pc*Pn*0.5*Pc*Pn*0.25*(1-(1-Py)^2)
  P5XOh<- (1-Pc)*0.5*Pc*(1-Pn)*0.5*Py + Pc*(1-Pn)*(1-Pc)*0.5*0.5*Py + 2*((1-Pc)*0.5*(1-Pc)*0.5*0.5*Py) +  (1-Pc)*0.5*Pc*Pn*0.25*Py +  Pc*Pn*0.5*(1-Pc)*0.25*Py
  P5XOw<- 0
  P5nv <-  Pc*Pn*0.5*Pc*(1-Pn)*0.5 +  Pc*Pn*0.5*(1-Pc)*0.25 +  Pc*(1-Pn)*Pc*Pn*0.25 +  (1-Pc)*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*Pc*(1-Pn)*0.5 +  Pc*Pn*0.5*(1-Pc)*0.25 +  Pc*(1-Pn)*Pc*Pn*0.25 +  (1-Pc)*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*Pc*Pn*0.25 +  (1-Pc)*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*(1-Pc)*0.25 +  (1-Pc)*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*(1-Pc)*0.25 +  Pc*Pn*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*Pc*Pn*0.25 +  Pc*Pn*0.5*Pc*Pn*0.25 
  
  #P6: XXh * XYH
  P6XXw<- 0
  P6XYw<- 0
  P6XXH<- (1-Pc)*0.25 + Pc*(1-Pn)*0.5 + Pc*Pn*0.25
  P6XXh<- (1-Pc)*0.25
  P6XYH<- (1-Pc)*0.25*(1-Py)^2 + Pc*(1-Pn)*0.5*(1-Py)^2 + Pc*Pn*0.25*(1-Py)^2
  P6XYh<- (1-Pc)*0.25*(1-Py)
  P6XOH<- (1-Pc)*0.25*(1-(1-Py)^2) + Pc*(1-Pn)*0.5*(1-(1-Py)^2) + Pc*Pn*0.25*(1-(1-Py)^2)
  P6XOh<- (1-Pc)*0.25*Py
  P6XOw<- 0
  P6nv <- 2*(Pc*Pn*0.25)
  
  #P7: XXH * XYw
  P7XXw<- 0
  P7XYw<- 0
  P7XXH<- 0
  P7XXh<- 0.5
  P7XYH<- 0
  P7XYh<- 0.5*(1-Py)
  P7XOH<- 0
  P7XOh<- 0.5*Py
  P7XOw<- 0
  P7nv <- 0
  
  #P8: XXH * XYh
  P8XXw<- 0
  P8XYw<- 0
  P8XXH<- Pc*(1-Pn)*0.5 + (1-Pc)*0.25 + Pc*Pn*0.25
  P8XXh<- (1-Pc)*0.25
  P8XYH<- Pc*(1-Pn)*0.5*(1-Py)^2 + (1-Pc)*0.25*(1-Py)^2 + Pc*Pn*0.25*(1-Py)^2
  P8XYh<- (1-Pc)*0.25*(1-Py)
  P8XOH<- Pc*(1-Pn)*0.5*(1-(1-Py)^2) + (1-Pc)*0.25*(1-(1-Py)^2) + Pc*Pn*0.25*(1-(1-Py)^2)
  P8XOh<- (1-Pc)*0.25*Py
  P8XOw<- 0
  P8nv <- 2*(Pc*Pn*0.25)
  
  #P9: XXH * XYH
  P9XXw<- 0
  P9XYw<- 0
  P9XXH<- 0.5
  P9XXh<- 0
  P9XYH<- 0.5*(1-Py)^2
  P9XYh<- 0
  P9XOH<- 0.5*(1-(1-Py)^2)
  P9XOh<- 0
  P9XOw<- 0
  P9nv <- 0
  
  #P10: XOh * XYw
  P10XXw<- (1-Pc)*Xpass*0.25
  P10XYw<- (1-Pc)*Xpass*0.25
  P10XXH<- 0
  P10XXh<- (1-Pc)*Xpass*0.25 + Pc*(1-Pn)*Xpass*0.5 + Pc*Pn*Xpass*0.25
  P10XYH<- 0
  P10XYh<- (1-Pc)*Xpass*0.25*(1-Py) + Pc*(1-Pn)*Xpass*0.5*(1-Py) + Pc*Pn*Xpass*0.25*(1-Py)
  P10XOH<- 0
  P10XOh<- (1-Pc)*Xpass*0.25*Py + Pc*(1-Pn)*Xpass*0.5*Py + (1-Pc)*(1-Xpass)*0.25 + Pc*(1-Pn)*(1-Xpass)*0.5 + Pc*Pn*Xpass*0.25*Py + Pc*Pn*(1-Xpass)*0.25
  P10XOw<- (1-Pc)*(1-Xpass)*0.25
  P10nv <- 2*(Pc*Pn*Xpass*0.25) + Pc*Pn*(1-Xpass)*0.25 + 2*((1-Pc)*(1-Xpass)*0.25) + Pc*(1-Pn)*(1-Xpass)*0.5 + Pc*Pn*(1-Xpass)*0.25 + Pc*Pn*(1-Xpass)*0.25
  
  #P11: XOh * XYh
  P11XXw<- (1-Pc)*Xpass*0.125*(1-Pc)
  P11XYw<- (1-Pc)*Xpass*0.125*(1-Pc)
  P11XXH<- (1-Pc)*Xpass*0.125*(1-Pc) + Pc*(1-Pn)*Xpass*(1-Pc)*0.25 + Pc*Pn*Xpass*0.125*(1-Pc) + Pc*(1-Pn)*Xpass*Pc*(1-Pn)*0.5 + Pc*(1-Pn)*Xpass*Pc*Pn*0.5*0.5 + Pc*Pn*Xpass*0.5*Pc*(1-Pn)*0.5 + Pc*Pn*Xpass*0.5*Pc*Pn*0.25 + (1-Pc)*Xpass*0.5*Pc*(1-Pn)*0.5 + (1-Pc)*Xpass*0.5*Pc*Pn*0.25
  P11XXh<- 2*((1-Pc)*Xpass*0.5*(1-Pc)*0.25) + Pc*(1-Pn)*Xpass*(1-Pc)*0.25 + Pc*Pn*Xpass*0.5*(1-Pc)*0.25 + (1-Pc)*Xpass*0.5*Pc*(1-Pn)*0.5 + (1-Pc)*Xpass*0.125*Pc*Pn
  P11XYH<- (1-Pc)*Xpass*0.125*(1-Pc)*(1-Py)^2 + Pc*(1-Pn)*Xpass*(!-Pc)*0.25*(1-Py)^2 + Pc*Pn*Xpass*0.125*(1-Pc)*(1-Py)^2 + Pc*(1-Pn)*Xpass*Pc*(1-Pn)*0.5*(1-Py)^2 + Pc*Pn*Xpass*0.25*Pc*(1-Pn)*(1-Py)^2 + Pc*(1-Pn)*Xpass*Pc*Pn*0.25*(1-Py)^2 + Pc*Pn*Xpass*0.125*Pc*Pn*(1-Py)^2 + (1-Pc)*Xpass*0.25*Pc*(1-Pn)*(1-Py)^2 + (1-Pc)*Xpass*0.125*Pc*Pn*(1-Py)^2
  P11XYh<- 2*((1-Pc)*Xpass*0.5*(1-Pc)*0.5*0.5*(1-Py)) + Pc*(1-Pn)*Xpass*(1-Pc)*0.25*(1-Py) + Pc*Pn*Xpass*0.125*(1-Pc)*(1-Py) + (1-Pc)*Xpass*0.25*Pc*(1-Pn)*(1-Py) + (1-Pc)*Xpass*0.125*Pc*Pn*(1-Py)
  P11XOH<- (1-Pc)*Xpass*0.125*(1-Pc)*(1-(1-Py)^2) + (1-Pc)*(1-Xpass)*0.125*(1-Pc) +  Pc*(1-Pn)*Xpass*(1-Pc)*0.25*(1-(1-Py)^2) + Pc*Pn*Xpass*0.125*(1-Pc)*(1-(1-Py)^2) + Pc*(1-Pn)*(1-Xpass)*(1-Pc)*0.25 + Pc*Pn*(1-Xpass)*0.125*(1-Pc) + (1-Pc)*Xpass*0.25*Pc*(1-Pn)*(1-(1-Py)^2) + (1-Pc)*Xpass*0.125*Pc*Pn*(1-(1-Py)^2) + (1-Pc)*(1-Xpass)*0.125*Pc*(1-Pn) + (1-Pc)*(1-Xpass)*0.125*Pc*Pn + Pc*(1-Pn)*Xpass*Pc*(1-Pn)*0.5*(1-(1-Py)^2) +  Pc*Pn*Xpass*0.25*Pc*(1-Pn)*(1-(1-Py)^2) + Pc*(1-Pn)*Xpass*Pc*Pn*0.25*(1-(1-Py)^2) + Pc*Pn*Xpass*0.125*Pc*Pn*(1-(1-Py)^2) +  Pc*(1-Pn)*(1-Xpass)*Pc*(1-Pn)*0.5 + Pc*Pn*(1-Xpass)*0.25*Pc*(1-Pn) +  Pc*(1-Pn)*(1-Xpass)*Pc*Pn*0.25 +  Pc*Pn*(1-Xpass)*0.125*Pc*Pn
  P11XOh<- 2*((1-Pc)*Xpass*0.125*(1-Pc)*Py) + 2*((1-Pc)*(1-Xpass)*0.125*(1-Pc)*Py) + Pc*(1-Pn)*Xpass*(1-Pc)*0.25*Py + Pc*Pn*Xpass*0.125*(1-Pc)*Py + Pc*(1-Pn)*(1-Xpass)*(1-Pc)*0.25 + Pc*Pn*(1-Xpass)*(1-Pc)*0.25 + (1-Pc)*Xpass*0.25*Pc*(1-Pn)*Py + (1-Pc)*Xpass*0.125*Pc*Pn*Py + (1-Pc)*(1-Xpass)*0.25*Pc*(1-Pn) + (1-Pc)*(1-Xpass)*0.25*Pc*Pn
  P11XOw<- (1-Pc)*(1-Xpass)*0.125*(1-Pc)
  P11nv <-   Pc*Pn*Xpass*0.125*(1-Pc) +   Pc*(1-Pn)*Xpass*Pc*Pn*0.5*0.5 +   Pc*Pn*Xpass*0.5*Pc*(1-Pn)*0.5 +   Pc*Pn*Xpass*0.5*Pc*Pn*0.25 +   (1-Pc)*Xpass*0.5*Pc*Pn*0.25 +   Pc*Pn*Xpass*0.125*(1-Pc) +   (1-Pc)*Xpass*0.125*Pc*Pn +   Pc*Pn*Xpass*0.125*(1-Pc) +   Pc*Pn*Xpass*0.25*Pc*(1-Pn) +   Pc*(1-Pn)*Xpass*Pc*Pn*0.25 +   Pc*Pn*Xpass*0.125*Pc*Pn +   (1-Pc)*Xpass*0.125*Pc*Pn +   Pc*Pn*Xpass*0.125*(1-Pc) +  (1-Pc)*Xpass*0.125*Pc*Pn +   Pc*Pn*Xpass*0.125*(1-Pc) +  Pc*Pn*(1-Xpass)*0.125*(1-Pc) +   (1-Pc)*Xpass*0.125*Pc*Pn  +   (1-Pc)*(1-Xpass)*0.125*Pc*Pn +   Pc*Pn*Xpass*0.25*Pc*(1-Pn) +    Pc*(1-Pn)*Xpass*Pc*Pn*0.25  +   Pc*Pn*Xpass*0.125*Pc*Pn +   Pc*Pn*(1-Xpass)*0.25*Pc*(1-Pn) +   Pc*(1-Pn)*(1-Xpass)*Pc*Pn*0.25 +   Pc*Pn*(1-Xpass)*0.125*Pc*Pn +   Pc*Pn*Xpass*0.125*(1-Pc) +   Pc*Pn*(1-Xpass)*(1-Pc)*0.25 +   (1-Pc)*Xpass*0.125*Pc*Pn +   (1-Pc)*(1-Xpass)*0.125*Pc*Pn +  (1-Pc)*(1-Xpass)*0.125*(1-Pc) +  Pc*(1-Pn)*(1-Xpass)*(1-Pc)*0.5*0.5 +    Pc*Pn*(1-Xpass)*0.125*(1-Pc) +    (1-Pc)*(1-Xpass)*0.25*Pc*(1-Pn) +    (1-Pc)*(1-Xpass)*0.125*Pc*Pn +    Pc*(1-Pn)*(1-Xpass)*Pc*(1-Pn)*0.5 + 2*((1-Pc)*(1-Xpass)*0.125*(1-Pc)) +   Pc*(1-Pn)*(1-Xpass)*(1-Pc)*0.25 +    Pc*Pn*(1-Xpass)*0.125*(1-Pc) +    (1-Pc)*(1-Xpass)*0.25*Pc*(1-Pn) +    (1-Pc)*(1-Xpass)*0.25*Pc*Pn +   (1-Pc)*(1-Xpass)*0.125*(1-Pc)  +Pc*Pn*(1-Xpass)*Pc*(1-Pn)*0.25 +  Pc*(1-Pn)*(1-Xpass)*Pc*Pn*0.25 + Pc*Pn*(1-Xpass)*Pc*Pn*0.125
  
  #P12: XOh * XYH
  P12XXw<- 0
  P12XYw<- 0
  P12XXH<- (1-Pc)*Xpass*0.25 + Pc*(1-Pn)*Xpass*0.5 + Pc*Pn*Xpass*0.25
  P12XXh<- (1-Pc)*Xpass*0.25
  P12XYH<- (1-Pc)*Xpass*0.25*(1-Py)^2 + Pc*(1-Pn)*Xpass*0.5*(1-Py)^2 + Pc*Pn*Xpass*0.25*(1-Py)^2
  P12XYh<- (1-Pc)*Xpass*0.25*(1-Py)
  P12XOH<- (1-Pc)*Xpass*0.25*(1-(1-Py)^2) + Pc*(1-Pn)*Xpass*0.5*(1-(1-Py)^2) + Pc*Pn*Xpass*0.25*(1-(1-Py)^2) + (1-Pc)*(1-Xpass)*0.25 + Pc*(1-Pn)*(1-Xpass)*0.5 + Pc*Pn*(1-Xpass)*0.25
  P12XOh<- (1-Pc)*Xpass*0.25*Py + (1-Pc)*0.25*(1-Xpass)
  P12XOw<- 0
  P12nv <- Pc*Pn*Xpass*0.25 + Pc*Pn*Xpass*0.25 + Pc*Pn*(1-Xpass)*0.25 + (1-Pc)*(1-Xpass)*0.25 + Pc*(1-Pn)*(1-Xpass)*0.5 + Pc*Pn*0.25*(1-Xpass) + (1-Pc)*0.25*(1-Xpass) +Pc*Pn*0.25*(1-Xpass)
  
  #P13: XOH * XYw
  P13XXw<- 0
  P13XYw<- 0
  P13XXH<- 0
  P13XXh<- Xpass*0.5
  P13XYH<- 0
  P13XYh<- Xpass*0.5*(1-Py)
  P13XOH<- 0
  P13XOh<- Xpass*0.5*Py + (1-Xpass)*0.5
  P13XOw<- 0
  P13nv <- (1-Xpass)*0.5
  
  #P14: XOH * XYh
  P14XXw<- 0
  P14XYw<- 0
  P14XXH<- Xpass*(1-Pc)*0.25 + Xpass*Pc*(1-Pn)*0.5 + Xpass*Pc*Pn*0.25
  P14XXh<- Xpass*(1-Pc)*0.25
  P14XYH<- Xpass*(1-Pc)*0.25*(1-Py)^2 + Xpass*Pc*(1-Pn)*0.5*(1-Py)^2 + Xpass*Pc*Pn*0.25*(1-Py)^2
  P14XYh<- Xpass*(1-Pc)*0.25*(1-Py)
  P14XOH<- Xpass*(1-Pc)*0.25*(1-(1-Py)^2) + Xpass*Pc*(1-Pn)*0.5*(1-(1-Py)^2) + Xpass*Pc*Pn*0.25*(1-(1-Py)^2) + (1-Xpass)*(1-Pc)*0.25 + (1-Xpass)*Pc*(1-Pn)*0.5 + (1-Xpass)*Pc*Pn*0.25
  P14XOh<- Xpass*(1-Pc)*0.25*Py + (1-Xpass)*(1-Pc)*0.25 
  P14XOw<- 0
  P14nv <- Xpass*Pc*Pn*0.25 + Xpass*Pc*Pn*0.25 + (1-Xpass)*Pc*Pn*0.25 + (1-Xpass)*(1-Pc)*0.25 + (1-Xpass)*Pc*(1-Pn)*0.5 + (1-Xpass)*Pc*Pn*0.25 + (1-Xpass)*(1-Pc)*0.25
  
  #P15: XOH * XYH
  P15XXw<- 0
  P15XYw<- 0
  P15XXH<- Xpass*0.5
  P15XXh<- 0
  P15XYH<- Xpass*0.5*(1-Py)^2
  P15XYh<- 0
  P15XOH<- (1-Xpass)*0.5 + (Xpass)*0.5*(1-(1-Py)^2)
  P15XOh<- 0
  P15XOw<- 0
  P15nv <- (1-Xpass)*0.5
  
  #P16: XOw * XYw
  P16XXw<- Xpass*0.5
  P16XYw<- Xpass*0.5
  P16XXH<- 0
  P16XXh<- 0
  P16XYH<- 0
  P16XYh<- 0
  P16XOH<- 0
  P16XOh<- 0
  P16XOw<- (1-Xpass)*0.5
  P16nv <- (1-Xpass)*0.5
  
  #P17: XOw * XYH
  P17XXw<- 0
  P17XYw<- 0
  P17XXH<- 0
  P17XXh<- Xpass*0.5
  P17XYH<- 0
  P17XYh<- Xpass*0.5*(1-Py)
  P17XOH<- 0
  P17XOh<- Xpass*0.5*Py + (1-Xpass)*0.5
  P17XOw<- 0
  P17nv <- (1-Xpass)*0.5
  
  #P18: XOw * XYh
  P18XXw<- Xpass*(1-Pc)*0.25
  P18XYw<- Xpass*(1-Pc)*0.25
  P18XXH<- 0
  P18XXh<- Xpass*(1-Pc)*0.25 + Xpass*Pc*(1-Pn)*0.5 + Xpass*Pc*Pn*0.25
  P18XYH<- 0
  P18XYh<- Xpass*(1-Pc)*0.25*(1-Py) + Xpass*Pc*(1-Pn)*0.5*(1-Py) + Pc*Pn*0.25*(1-Py)
  P18XOH<- 0
  P18XOh<- Xpass*(1-Pc)*0.25*Py + Xpass*Pc*(1-Pn)*0.5*Py + Xpass*Pc*Pn*0.25*Py + (1-Xpass)*(1-Pc)*0.25 + (1-Xpass)*Pc*(1-Pn)*0.5 + (1-Xpass)*Pc*Pn*0.25
  P18XOw<- (1-Xpass)*(1-Pc)*0.25
  P18nv <- Xpass*Pc*Pn*0.25 + Pc*Pn*0.25 + (1-Xpass)*Pc*Pn*0.25 + (1-Xpass)*(1-Pc)*0.25 + (1-Xpass)*Pc*(1-Pn)*0.5 + (1-Xpass)*Pc*Pn*0.25 + (1-Xpass)*(1-Pc)*0.25
  
  ##############################################################################
  temp<-(ode(func = Yshred, y = yn, times = 0:(weeks-1), parms = tparms[,i], method = "iteration"))
  temp.df<-as.data.frame(temp[1:weeks,(ncol(temp)-7):ncol(temp)])
  #write.csv(temp.df, paste0(outloc,"/",outfolder,"_",i,".csv"), sep=",", row.names = F)
  ########################
  #Compile results into shared FBM matrix
  bir0<-which(temp.df[,1]==0)  #pull out time steps with no births
  if(length(bir0)>2){    #if there are 3 or more times than check for consecutives
    for(k in 1:(length(bir0)-2)){
      if(bir0[k+1]==bir0[k]+1 & bir0[k+2]==bir0[k]+2){   #if there are 3 consecutive times with 0 births
        runResults[i,]<-c(tparms[22,i],i, 1, bir0[k], (1-(temp.df[weeks,5]/temp.df[23,5])))
        {break}
      } else {
        runResults[i,]<-c(tparms[22,i],i,0,0,(1-(temp.df[weeks,5]/temp.df[23,5])))
      }
    }
  } else{
    runResults[i,]<-c(tparms[22,i],i,0,0,(1-(temp.df[weeks,5]/temp.df[23,5])))
  }
  rm(temp)
  
  if(i%%10000==0){
    big_write(runResults, paste0(outloc,"/",outfolder,"_progRes_",i,".csv"), every_nrow = 1)
  }
}

#Stop cluster
parallel::stopCluster(cl = my.cluster)
#Record end time
end.time <- Sys.time()  
time.taken <- end.time - start.time
print(time.taken)

big_write(runResults, paste0(outloc,"/",outfolder,"_runResults.csv"), every_nrow = 1)
