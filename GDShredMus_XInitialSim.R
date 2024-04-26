### Title: Script to run initial simulations for X-Shredder gene drive
### Author: Matthew Combs
### Date: 26APR2024
### Notes:
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
outfolder<-"XShred_InitialSim"
basepath.out<-"~/outfiles/"
outloc<-paste0(basepath.out, outfolder)
dir.create(outloc)


##############################################################################
#Set lattice size, time steps, iterations
NL<-31 #Length/Width of square lattice of patches
weeks<- 742 #30 years + 22 timesteps (burn-in)
iterations<-500

##############################################################################
#Create parameter matrix
eFactorial<- expand.grid(Factor1 = seq(0.5, 0.5, length.out=1),   #Development rate
                          Factor2 = seq(0.07, 0.07, length.out=1),     #Mortality rate
                          Factor3 = seq(6,6,length.out=1),             #Litter size
                          Factor4 = seq(0.44, 0.44, length.out=1),     #Female mating rate
                          Factor5 = seq(50, 50, length.out=1),         #Patch carrying capacity
                          Factor6 = seq(1, 4, length.out=4),           #MaxAndry: max male mates per female
                          Factor7 = seq(1, 4, length.out=4),           #MaxGyny: max female mates per male
                          Factor8 = seq(2, 2, length.out=1),           #Dispersal kernel shape
                          Factor9 = seq(3.0, 0.6, length.out=7),       #Dispersal kernel rate (influences avg distance)
                          Factor10 = c(0.0125, 0.025, 0.05, 0.1, 0.15, 0.2), #Juvenile dispersal probability
                          Factor11 = seq(0.0, 0.0, length.out=1),      #Adult dispersal probability
                          Factor12 = seq(1000, 1000, length.out=1),    #Intro size
                          Factor13 = seq(11, 11, length.out=1),        #GD Supplementation schedule
                          Factor14 = seq(2, 2, length.out=1),          #Placement on patch lattice
                          Factor15 = seq(1.0, 1.0, length.out=1),      #Drive efficiency/Inverse resistance rate
                          Factor16 = seq(0.75, 0.75, length.out=1),   #Sperm Competitiveness
                          Factor17 = seq(2, 2, length.out=1)           #Polyandry Rate
)

vec<-c(1:672)
eFactorial[1:nrow(eFactorial),(ncol(eFactorial)+1)]<-vec  
eFactorial <- as.data.frame(t(eFactorial))
tparms<- matrix(data = apply(eFactorial, 2, function(x) rep(x, iterations)), ncol = ncol(eFactorial)*iterations)

#Setup matrix to write basic results
runResults<-FBM(ncol(tparms),5)

##############################################################################
#Additional run setup

#Null mating matrix, used in randBirthXS function
b0 <-rep.int(0, NL*NL*3)
Birthnull<-matrix(b0,ncol=(3),nrow=NL*NL,byrow = T) 

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
  ########################
  #Starting pops based on carrying capacity
  #Note: placement of GDMO in center cell only works if NL is odd number
  #Males
  XY0 <-matrix((floor(tparms[5,i] * 0.5)),NL,NL)
  XD0 <-matrix(0,NL,NL)
  #Females
  XX0 <-matrix((floor(tparms[5,i] * 0.5)),NL,NL)
  #Juv Males
  jXY0 <-matrix(0,NL,NL)
  jXD0 <-matrix(0,NL,NL)
  #Juv Females
  jXX0 <-matrix(0,NL,NL)
  
  #Calculate summary values
  patchTotal <- XY0+XD0+XX0+jXY0+jXD0+jXX0
  Total0 <- sum(patchTotal)
  GDMales <- XD0/(XD0+XY0)
  #Empty set of 0s for births at first time step
  birthsempty<-matrix(0,NL,NL)
  #Average percentage of indivduals that are GDMO across patches
  Out10 <- sum(XD0)/sum(Total0) 
  #Overall abundance of all individuals
  Out30 <- sum(Total0)
  #Sex Ratio
  sexRatio<- (XY0+XD0+jXY0+jXD0) / Total0
  #Inv Frq
  Present <- ifelse(XD0<1,0,1)
  Present <- as.vector(Present)
  InvFrq <- sum(Present)/(NL^2)
  #Extinction Frq
  Extinct <- ifelse(patchTotal<1, 1, 0)
  ExtFrq <- sum(Extinct)/(NL^2)
  #Geno Frqs
  fqnXD0<-sum(XD0+jXD0)/Total0
  fqnXY0<-sum(XY0+jXY0)/Total0
  fqnXX0<-sum(XX0+jXX0)/Total0
  
  #Combine into list of matrices to initialize model
  yn <-c(XY0,XD0,XX0,jXY0,jXD0,jXX0,
         patchTotal, sexRatio, birthsempty, GDMales,
         0,0,Out10, Out30,
         InvFrq,ExtFrq,fqnXD0,fqnXY0,fqnXX0,1)
  
  #Establish mating system parameters for specific iteration
  PolyandrySeq<-seq(1,tparms[6,i],1)   #setup for polyandry 
  AndryRate<-tparms[17,i]
  if(AndryRate==1){
    PolyandryProbs<-PolyandryProbsList[[(tparms[6,i])]]
  } else if(AndryRate==2){
    PolyandryProbs<-PolyandryProbsList[[ 4 + (tparms[6,i]) ]]
  } else if(AndryRate==3){
    PolyandryProbs<-PolyandryProbsList[[ 8 + (tparms[6,i]) ]]
  }
  #Genotype calcs based on Pc param, used for reproduction
  Pc<-tparms[15,i]
  P1XX<-0.5
  P1XY<-0.5
  P2XX<-(1-Pc)*0.5
  P2XD<-Pc+(1-Pc)*0.5
  ########################
  temp<-(ode(func = Xshred.v2, y = yn, times = 0:(weeks-1), parms = tparms[,i], method = "iteration"))
  temp.df<-as.data.frame(temp[1:weeks,(ncol(temp)-9):ncol(temp)])
  #write.csv(temp.df, paste0(outloc,"/",outfolder,"_",i,".csv"), sep=",", row.names = F)
 
   ########################
  #Compile results into shared FBM matrix
  bir0<-which(temp.df[,1]==0)  #pull out time steps with no births
  if(length(bir0)>2){    #if there are 3 or more times than check for consecutives
    for(k in 1:(length(bir0)-2)){
      if(bir0[k+1]==bir0[k]+1 & bir0[k+2]==bir0[k]+2){   #if there are 3 consecutive times with 0 births
        runResults[i,]<-c(tparms[18,i],i, 1, bir0[k], (1-(temp.df[weeks,4]/temp.df[23,4])))
        {break}
      } else {
        runResults[i,]<-c(tparms[18,i],i,0,0,(1-(temp.df[weeks,4]/temp.df[23,4])))
      }
    }
  } else{
    runResults[i,]<-c(tparms[18,i],i,0,0,(1-(temp.df[weeks,4]/temp.df[23,4])))
  }
  rm(temp)
  
  if(i%%5000==0){
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
