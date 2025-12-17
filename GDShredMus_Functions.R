### Title: Functions for simulating X-Shredder and Y-Shredder gene drives for house mouse metapopulations
### Author: Matthew Combs
### Date: December 2025
### Notes:

#X-shredder gene drive dynamics
Xshred <- function(t, y, parms, signal){
  
  ################################################################################  
  #Setup: Load Parameters and initital abundances
  N<-NL
  g=parms[[1]]            #Development rate
  u=parms[[2]]            #Mortality rate
  k=parms[[3]]            #Litter size
  v=parms[[4]]            #Female mating rate
  CPop=parms[[5]]         #Carrying capacity of each patch
  Andry=parms[[6]]        #Max male mates per female
  Gyny=parms[[7]]         #Max female mates per male
  kerShape=parms[[8]]     #Dispersal kernel shape
  kerRate=parms[[9]]      #dispersal kernel rate
  dJ=parms[[10]]          #Proportion of juveniles dispersing
  dA=parms[[11]]          #Proportion of adults dispersing
  gdInt=parms[[12]]       #number of gene drive individuals introduced
  gdSupp=parms[[13]]      #release schedule (every 1-10yrs), or once
  intLoc=parms[[14]]      #release all in middle (1) or each ind randomly (2)
  e=parms[[15]]           #prob of HEG1 cutting at homologous locus
  spD=parms[[16]]         #sperm competition disadvantage for drive males (XD)
  Dyr=parms[[18]]         #how many iterations between dispersal events
  
  #Setup: Create vector of N from vector loaded into model
  XY <- y[(0*(N^2) + 1):(1*(N^2))] 
  XD <- y[(1*(N^2) + 1):(2*(N^2))]
  XX <- y[(2*(N^2) + 1):(3*(N^2))]
  jXY <- y[(3*(N^2) + 1):(4*(N^2))]
  jXD <- y[(4*(N^2) + 1):(5*(N^2))]
  jXX <- y[(5*(N^2) + 1):(6*(N^2))]
  time.step<-y[10*(N^2)+10]
  
  ################################################################################  
  #GD Release at specific time steps
  if(gdSupp==11){
    if(time.step==7){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- gdInt #Add specific n of GD ind to central patch 
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
  }
  
  if(gdSupp<11) {
    if(time.step==7){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- gdInt #Always add at the beginning
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
    
    if(time.step>7 & (((time.step-7)/6)%%gdSupp)==0){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- (XD[(N^2/2)+0.5]+gdInt) #Add new GD individuals to central patch
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
  }
  
  ################################################################################
  #Step 1: Mortality and Maturation
  
  #Adult mortality
  D1 <- (rbinom(n = XY, size=XY, prob = u))
  D2 <- (rbinom(n = XD, size=XD, prob = u))
  D3 <- (rbinom(n = XX, size=XX, prob = u))
  
  #Juvenile mortality
  D4 <- (rbinom(n = jXY, size=jXY, prob = u))
  D5 <- (rbinom(n = jXD, size=jXD, prob = u))
  D6 <- (rbinom(n = jXX, size=jXX, prob = u))
  
  #Update juvenile population after mortality
  jXY <- jXY - D4
  jXD <- jXD - D5
  jXX <- jXX - D6
  jXY <- ifelse(jXY<0,0,jXY) 
  jXD <- ifelse(jXD<0,0,jXD) 
  jXX <- ifelse(jXX<0,0,jXX) 
  jXY[is.na(jXY)] = 0
  jXD[is.na(jXD)] = 0
  jXX[is.na(jXX)] = 0
  
  #Juvenile maturation to adult
  J1 <- (rbinom(n = jXY, size=jXY, prob = g))
  J2 <- (rbinom(n = jXD, size=jXD, prob = g))
  J3 <- (rbinom(n = jXX, size=jXX, prob = g))
  
  #Update adult subpopulations
  mXY <- XY + J1 - D1
  mXD <- XD + J2 - D2
  mXX <- XX + J3 - D3
  mXY <- ifelse(mXY<0,0,mXY) 
  mXD <- ifelse(mXD<0,0,mXD) 
  mXX <- ifelse(mXX<0,0,mXX) 
  mXY[is.na(mXY)] = 0
  mXD[is.na(mXD)] = 0
  mXX[is.na(mXX)] = 0
  
  #Update juvenile population after maturation
  jXY <- jXY - J1
  jXD <- jXD - J2
  jXX <- jXX - J3
  jXY <- ifelse(jXY<0,0,jXY) 
  jXD <- ifelse(jXD<0,0,jXD) 
  jXX <- ifelse(jXX<0,0,jXX) 
  jXY[is.na(jXY)] = 0
  jXD[is.na(jXD)] = 0
  jXX[is.na(jXX)] = 0
  
  ################################################################################
  #Step 2: Adult mating
  
  #A. Mating population established (mating males "mM" and mating females "mF")
  mM <- mXY + mXD   #Male mating population
  mF <- v*(mXX) #Female Mating Population 
  
  #Number mated females given max number of mates per males and females
  #Requires polyandry probabilities (likelihood of n mates for females for 1:max)
  
  if(Andry==1){
    FnMated<-mapply(minCheck, round(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq)), mF)
    FnMated[is.na(FnMated)] = 0 
  }else{
    FnMated<-mapply(minCheck, round(colSums(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq))), mF)
    FnMated[is.na(FnMated)] = 0 
  }
  
  #B. Total number of births in each patch w density dependence, given avg litter size "k"
  sumLit<-(randpoisLitter(FnMated, k, FnMated))  #assumes min=2, max=9
  B<- sumLit*(1-((mXY + mXD + mXX)/CPop))
  B<- ifelse(B<0, 0, B)
  B <- (randpois(dsubpop = B, lambda = 1, subpop = B))
  
  #C. Mating fractions among all adult genotype pairs, assuming random mating within spatial patches
  allmate<- (mXX * (mXY + mXD))
  P1<- (mXX * mXY) /allmate
  P2<- (mXX * mXD) /allmate
  
  #D. calculate probability of offspring genotypes
  oXY<- P1*P1XY
  oXD<- P2*P2XD
  oXX<- P1*P1XX + P2*P2XX
  
  #E. Adjust offspring probabilities due to fitness costs. GD inds have reduced likelihood of birth
  if(length(PolyandrySeq)>1){
    oXD <- oXD-(oXD*(1-spD)*(sum(PolyandryProbs[-1])))
  }
  
  #F. Establish matrix of mating probabilities for randBirth2 function
  pBirth <-matrix(c(oXY,oXD,oXX),nrow=N*N,ncol=3,byrow=FALSE) #Male Bias. 
  
  #G. Establish final integer number of offpsring of each genotype
  Births <- randBirthXS(Birth = Birthnull, pBirth = pBirth, B = B)
  B1 <- Births[,1]
  B2 <- Births[,2]
  B3 <- Births[,3]
  
  
  ################################################################################
  #Step 3: Adult dispersal

  if((time.step-8)%%Dyr==0){
    
    M1 <- (rbinom(n = mXY, size=mXY, prob = dA))
    M2 <- (rbinom(n = mXD, size=mXD, prob = dA))
    M3 <- (rbinom(n = mXX, size=mXX, prob = dA))
    
    #Returns vector of new immigrants in each patch
    #Needs shape and rate parameters
    #Needs pSize parameter for informing range of potential destination patches
    tM1 <- (natMigration2(nMig=M1, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM2 <- (natMigration2(nMig=M2, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM3 <- (natMigration2(nMig=M3, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    
    #Update adult pops
    nXY <- mXY + tM1 - M1
    nXD <- mXD + tM2 - M2
    nXX <- mXX + tM3 - M3
    nXY <- ifelse(nXY<1, 0, nXY) #Ensure population does not go below 1. At least one needed to mate
    nXD <- ifelse(nXD<1, 0, nXD) 
    nXX <- ifelse(nXX<1, 0, nXX) 
  } else{
    nXY <- mXY
    nXD <- mXD
    nXX <- mXX
    nXY <- ifelse(nXY<1, 0, nXY) #Ensure population does not go below 1. At least one needed to mate
    nXD <- ifelse(nXD<1, 0, nXD) 
    nXX <- ifelse(nXX<1, 0, nXX) 
    
  }
  
  ##################################################################################
  #Step 4: Juvenile Dispersal
  
  if((time.step-8)%%Dyr==0){
    
    #Calc number of dispersers
    M4 <- (rbinom(n = B1, size=B1, prob = dJ))
    M5 <- (rbinom(n = B2, size=B2, prob = dJ))
    M6 <- (rbinom(n = B3, size=B3, prob = dJ))
    
    #Return number of new immigrants in each patch
    tM4 <- (natMigration2(nMig=M4, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM5 <- (natMigration2(nMig=M5, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM6 <- (natMigration2(nMig=M6, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    
    #Update Juvenile populations from births and dispersal
    njXY <- jXY + B1 - M4 + tM4 
    njXD <- jXD + B2 - M5 + tM5
    njXX <- jXX + B3 - M6 + tM6
    
    njXY <- ifelse(njXY<0, 0,njXY)
    njXD <- ifelse(njXD<0, 0,njXD)
    njXX <- ifelse(njXX<0, 0,njXX)
    njXY[is.na(njXY)] = 0
    njXD[is.na(njXD)] = 0
    njXX[is.na(njXX)] = 0
  } else {
    njXY <- jXY + B1 
    njXD <- jXD + B2
    njXX <- jXX + B3
    
    njXY <- ifelse(njXY<0, 0,njXY)
    njXD <- ifelse(njXD<0, 0,njXD)
    njXX <- ifelse(njXX<0, 0,njXX)
    njXY[is.na(njXY)] = 0
    njXD[is.na(njXD)] = 0
    njXX[is.na(njXX)] = 0
  }
  ################################################################################
  #Step 5: Record important results
  
  #####
  #Sum births
  totBirth<-sum(B)
  
  #Step 5: Calculate end of time-step totals
  patchTotal <- nXY+nXD+nXX+njXY+njXD+njXX
  PopTot <- sum(patchTotal) #Global population
  
  GDMO <- nXD  #not including nRD because dead end individuals that wont drive themselves
  GDMOTot <- sum(GDMO)   #Global GDMO
  
  #Binary vector indicating if GDMOs are present (1) or not (0) in each patch
  Present <- ifelse(GDMO<1, 0,1)
  Present <- as.vector(Present)
  
  #Freq of gd allele in male pop
  GDMales<- (nXD)/(nXY+nXD)
  
  #genotype total for animation
  fqnXD<-sum(nXD+njXD)/PopTot
  fqnXY<-sum(nXY+njXY)/PopTot
  fqnXX<-sum(nXX+njXX)/PopTot
  
  #Sex Ratio in each patch (>0.5 is skewed male, <0.5 is skewed female)
  sexRatio<- (nXY+nXD+njXY+njXD) / patchTotal
  
  #Number of births in each patch
  #already contained in vector B
  
  #Out1: global frequency of gdmo individuals in context of total population
  #      vector of size NL*NL
  GDFrq <- GDMOTot/PopTot #Global Frequency
  GDFrq[is.na(GDFrq)] = 0
  
  #Out3: frequency of cells with at least 1 gdmo adult present (invaded)
  InvFrq <- sum(Present)/(N*N)
  
  #Out4: frequency of patches that have 0 individuals wt or gdmo (extinct)
  Extinct <- ifelse(patchTotal<1, 1, 0)
  ExtFrq <- sum(Extinct)/(N*N)
  
  #Tracking the time step...maybe not necessary
  ntime.step <-(time.step+1)
  
  dx <- c(nXY,nXD,nXX,njXY,njXD,njXX,
          patchTotal, sexRatio, B, GDMales,
          totBirth, GDMOTot, GDFrq, PopTot,
          InvFrq, ExtFrq, fqnXD, fqnXY, fqnXX, ntime.step) #vector of populations.
  
  #######################################
  ##combine results into a single vector dx
  ######################################
  list(dx)  
}

#Y-Shredder gene drive dynamics
Yshred <- function(t, y, parms, signal){
  
  ################################################################################  
  #Setup: Load Parameters and initital abundances
  
  N<-NL
  g=parms[[1]]            #Development rate
  u=parms[[2]]            #Mortality rate
  k=parms[[3]]            #Litter size
  v=parms[[4]]            #Female mating rate
  CPop=parms[[5]]         #Carrying capacity of each patch
  Andry=parms[[6]]        #Max male mates per female
  Gyny=parms[[7]]         #Max female mates per male
  kerShape=parms[[8]]     #Dispersal kernel shape
  kerRate=parms[[9]]      #dispersal kernel rate
  dJ=parms[[10]]          #Proportion of juveniles dispersing
  dA=parms[[11]]          #Proportion of adults dispersing
  gdInt=parms[[12]]       #number of gene drive individuals introduced
  gdSupp=parms[[13]]      #release schedule (every 1-10yrs), or once
  intLoc=parms[[14]]      #release all in middle (1) or each ind randomly (2)
  Pc=parms[[15]]          #prob of HEG1 cutting at homologous locus
  Pn=parms[[16]]          #prob of NHEJ after HEG1 cut (i.e. successful conversion)
  Py=parms[[17]]          #prob of HEG2 cutting at Y chromosome
  Xpass=parms[[18]]       #prob of XO female passing X chrom
  s=parms[[19]]           #Fitness cost for homozygous GDMO
  sh=parms[[20]]          #Fitness cost for heterozygous GDMO
  
  #Setup: Create vector of N from vector loaded into model
  XXw <- y[(0*(N^2) + 1): (1*(N^2))] 
  XYw <- y[(1*(N^2) + 1): (2*(N^2))]
  XOw <- y[(2*(N^2) + 1): (3*(N^2))]
  XXH <- y[(3*(N^2) + 1): (4*(N^2))]
  XXh <- y[(4*(N^2) + 1): (5*(N^2))]
  XYH <- y[(5*(N^2) + 1): (6*(N^2))]
  XYh <- y[(6*(N^2) + 1): (7*(N^2))]
  XOH <- y[(7*(N^2) + 1): (8*(N^2))]
  XOh <- y[(8*(N^2) + 1): (9*(N^2))]
  jXXw <- y[(9*(N^2) + 1): (10*(N^2))]
  jXYw <- y[(10*(N^2) + 1): (11*(N^2))]
  jXOw <- y[(11*(N^2) + 1): (12*(N^2))]
  jXXH <- y[(12*(N^2) + 1): (13*(N^2))]
  jXXh <- y[(13*(N^2) + 1): (14*(N^2))]
  jXYH <- y[(14*(N^2) + 1): (15*(N^2))]
  jXYh <- y[(15*(N^2) + 1): (16*(N^2))]
  jXOH <- y[(16*(N^2) + 1): (17*(N^2))]
  jXOh <- y[(17*(N^2) + 1): (18*(N^2))]
  time.step<-y[((22*(N^2))+8)]
  
  ################################################################################  
  #GD Release at specific time steps
  
  #Scheduled release of GD individuals (heterozygous females), once here, successively below
  if(gdSupp==11){
    if(time.step==7){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- gdInt #Add specific n of GD ind to central patch 
      }
      if(intLoc==2){              #Random placement
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
  }
  #First release when supplemental releases scheduled
  if(gdSupp<11) {
    if(time.step==7){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- gdInt #Always add at the beginning
      }
      if(intLoc==2){              #Random placement
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
    #Supplemental release every X years
    if(time.step>7 & (((time.step-7)/6)%%gdSupp)==0){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- (XXH[(N^2/2)+0.5]+gdInt) #Add new GD individuals to central patch
      }
      if(intLoc==2){
        for(i in 1:gdInt){              #Random placement
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
  }
  
  ################################################################################
  #Step 1: Mortality and Maturation
  
  #Determine number of adult mortalities
  D1 <- (rbinom(n = XXw, prob = u, size = XXw))
  D2 <- (rbinom(n = XYw, prob = u, size = XYw))
  D3 <- (rbinom(n = XXH, prob = u, size = XXH))
  D4 <- (rbinom(n = XXh, prob = u, size = XXh))
  D5 <- (rbinom(n = XYH, prob = u, size = XYH))
  D6 <- (rbinom(n = XYh, prob = u, size = XYh))
  D7 <- (rbinom(n = XOH, prob = u, size = XOH))
  D8 <- (rbinom(n = XOh, prob = u, size = XOh))
  D9 <- (rbinom(n = XOw, prob = u, size = XOw))
  
  #Determine number of juvenile mortalities
  D10 <- (rbinom(n = jXXw, prob = u, size = jXXw))
  D11 <- (rbinom(n = jXYw, prob = u, size = jXYw))
  D12 <- (rbinom(n = jXXH, prob = u, size = jXXH))
  D13 <- (rbinom(n = jXXh, prob = u, size = jXXh))
  D14 <- (rbinom(n = jXYH, prob = u, size = jXYH))
  D15 <- (rbinom(n = jXYh, prob = u, size = jXYh))
  D16 <- (rbinom(n = jXOH, prob = u, size = jXOH))
  D17 <- (rbinom(n = jXOh, prob = u, size = jXOh))
  D18 <- (rbinom(n = jXOw, prob = u, size = jXOw))
  
  #Update juvenile population after mortality
  jXXw <- jXXw - D10
  jXYw <- jXYw - D11
  jXXH <- jXXH - D12
  jXXh <- jXXh - D13
  jXYH <- jXYH - D14
  jXYh <- jXYh - D15
  jXOH <- jXOH - D16
  jXOh <- jXOh - D17
  jXOw <- jXOw - D18
  jXXw <- ifelse(jXXw<0,0,jXXw) 
  jXYw <- ifelse(jXYw<0,0,jXYw) 
  jXXH <- ifelse(jXXH<0,0,jXXH) 
  jXXh <- ifelse(jXXh<0,0,jXXh) 
  jXYH <- ifelse(jXYH<0,0,jXYH) 
  jXYh <- ifelse(jXYh<0,0,jXYh) 
  jXOH <- ifelse(jXOH<0,0,jXOH) 
  jXOh <- ifelse(jXOh<0,0,jXOh) 
  jXOw <- ifelse(jXOw<0,0,jXOw) 
  jXXw[is.na(jXXw)] = 0
  jXYw[is.na(jXYw)] = 0
  jXXH[is.na(jXXH)] = 0
  jXXh[is.na(jXXh)] = 0
  jXYH[is.na(jXYH)] = 0
  jXYh[is.na(jXYh)] = 0
  jXOH[is.na(jXOH)] = 0
  jXOh[is.na(jXOh)] = 0
  jXOw[is.na(jXOw)] = 0
  
  #Poisson of maturation
  J1 <- (rbinom(n = jXXw, prob = g, size = jXXw))
  J2 <- (rbinom(n = jXYw, prob = g, size = jXYw))
  J3 <- (rbinom(n = jXXH, prob = g, size = jXXH))
  J4 <- (rbinom(n = jXXh, prob = g, size = jXXh))
  J5 <- (rbinom(n = jXYH, prob = g, size = jXYH))
  J6 <- (rbinom(n = jXYh, prob = g, size = jXYh))
  J7 <- (rbinom(n = jXOH, prob = g, size = jXOH))
  J8 <- (rbinom(n = jXOh, prob = g, size = jXOh))
  J9 <- (rbinom(n = jXOw, prob = g, size = jXOw))
  
  #Update adult subpopulations
  mXXw <- XXw + J1 - D1
  mXYw <- XYw + J2 - D2
  mXXH <- XXH + J3 - D3
  mXXh <- XXh + J4 - D4
  mXYH <- XYH + J5 - D5
  mXYh <- XYh + J6 - D6
  mXOH <- XOH + J7 - D7
  mXOh <- XOh + J8 - D8
  mXOw <- XOw + J9 - D9
  mXXw <- ifelse(mXXw<0,0,mXXw) 
  mXYw <- ifelse(mXYw<0,0,mXYw) 
  mXXH <- ifelse(mXXH<0,0,mXXH) 
  mXXh <- ifelse(mXXh<0,0,mXXh) 
  mXYH <- ifelse(mXYH<0,0,mXYH) 
  mXYh <- ifelse(mXYh<0,0,mXYh)
  mXOH <- ifelse(mXOH<0,0,mXOH)
  mXOh <- ifelse(mXOh<0,0,mXOh)
  mXOw <- ifelse(mXOw<0,0,mXOw)
  mXXw[is.na(mXXw)] = 0
  mXYw[is.na(mXYw)] = 0
  mXXH[is.na(mXXH)] = 0
  mXXh[is.na(mXXh)] = 0
  mXYH[is.na(mXYH)] = 0
  mXYh[is.na(mXYh)] = 0  
  mXOH[is.na(mXOH)] = 0 
  mXOh[is.na(mXOh)] = 0 
  mXOw[is.na(mXOw)] = 0 
  
  #Update juvenile population after maturation
  jXXw <- jXXw - J1
  jXYw <- jXYw - J2
  jXXH <- jXXH - J3
  jXXh <- jXXh - J4
  jXYH <- jXYH - J5
  jXYh <- jXYh - J6
  jXOH <- jXOH - J7
  jXOh <- jXOh - J8
  jXOw <- jXOw - J9
  jXXw <- ifelse(jXXw<0,0,jXXw) 
  jXYw <- ifelse(jXYw<0,0,jXYw) 
  jXXH <- ifelse(jXXH<0,0,jXXH) 
  jXXh <- ifelse(jXXh<0,0,jXXh) 
  jXYH <- ifelse(jXYH<0,0,jXYH) 
  jXYh <- ifelse(jXYh<0,0,jXYh) 
  jXOH <- ifelse(jXOH<0,0,jXOH) 
  jXOh <- ifelse(jXOh<0,0,jXOh) 
  jXOw <- ifelse(jXOw<0,0,jXOw) 
  jXXw[is.na(jXXw)] = 0
  jXYw[is.na(jXYw)] = 0
  jXXH[is.na(jXXH)] = 0
  jXXh[is.na(jXXh)] = 0
  jXYH[is.na(jXYH)] = 0
  jXYh[is.na(jXYh)] = 0
  jXOH[is.na(jXOH)] = 0
  jXOh[is.na(jXOh)] = 0
  jXOw[is.na(jXOw)] = 0
  
  ################################################################################
  #Step 2: Adult mating
  
  #A. Mating population established (mating males "mM" and mating females "mF")
  mM <- mXYw + mXYH + mXYh                    #Male mating population
  mF <- v*(mXXw + mXXH + mXXh + mXOH + mXOh + mXOw)  #Female Mating Population 
  
  # Determine number mated females given max number of mates per males and females
  #Requires polyandry probabilities (likelihood of n mates for females for 1:max)
  if(Andry==1){
    FnMated<-mapply(minCheck, round(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq)), mF)
    FnMated[is.na(FnMated)] = 0 
  }else{
    FnMated<-mapply(minCheck, round(colSums(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq))), mF)
    FnMated[is.na(FnMated)] = 0 
  }
  
  #B. Total number of births in each patch w density dependence, given avg litter size "k"
  sumLit<-(randpoisLitter(FnMated, k, FnMated))  #assumes min=2, max=9
  B<- sumLit*(1-((mXXw + mXYw + mXXH + mXXh + mXYH + mXYh + mXOH + mXOh + mXOw)/CPop))
  B<- ifelse(B<0, 0,B)
  B <- (randpois(dsubpop = B, lambda = 1, subpop = B))

  #C. Mating fractions among all adult genotype pairs, assuming random mating within spatial patches
  allmate<- (mXXw * (mXYw + mXYH +mXYh)) + (mXXH * (mXYw + mXYH +mXYh)) + (mXXh * (mXYw + mXYH +mXYh)) + (mXOH * (mXYw + mXYH +mXYh)) + (mXOh * (mXYw + mXYH +mXYh)) + (mXOw * (mXYw + mXYH +mXYh))
  P1<-(mXXw * mXYw) / allmate
  P2<-(mXXw * mXYh) / allmate
  P3<-(mXXw * mXYH) / allmate
  P4<-(mXXh * mXYw) / allmate
  P5<-(mXXh * mXYh) / allmate
  P6<-(mXXh * mXYH) / allmate
  P7<-(mXXH * mXYw) / allmate
  P8<-(mXXH * mXYh) / allmate
  P9<-(mXXH * mXYH) / allmate
  P10<-(mXOh * mXYw) / allmate
  P11<-(mXOh * mXYh) / allmate
  P12<-(mXOh * mXYH) / allmate
  P13<-(mXOH * mXYw) / allmate
  P14<-(mXOH * mXYh) / allmate
  P15<-(mXOH * mXYH) / allmate
  P16<-(mXOw * mXYw) / allmate
  P17<-(mXOw * mXYH) / allmate
  P18<-(mXOw * mXYh) / allmate
  
  #D. Offspring probs sum of mating fraction and inheritance calculation for all possible parent combos that could produce offspring geno
  oXXw<- P1*P1XXw + P2*P2XXw + P4*P4XXw + P5*P5XXw + P10*P10XXw + P11*P11XXw + P15*P15XXw + P18*P18XXw
  oXYw<- P1*P1XYw + P2*P2XYw + P4*P4XYw + P5*P5XYw + P10*P10XYw + P11*P11XYw + P15*P15XYw + P18*P18XYw
  oXXH<- P2*P2XXH + P3*P3XXH + P4*P4XXH + P5*P5XXH + P6*P6XXH + P7*P7XXH + P8*P8XXH + P9*P9XXH + P10*P10XXH + P11*P11XXH + P12*P12XXH + P13*P13XXH + P14*P14XXH + P15*P15XXH + P17*P17XXH + P18*P18XXH
  oXXh<- P2*P2XXh + P3*P3XXh + P4*P4XXh + P5*P5XXh + P6*P6XXh + P7*P7XXh + P8*P8XXh + P10*P10XXh + P11*P11XXh + P12*P12XXh + P13*P13XXh + P14*P14XXh + P17*P17XXh + P18*P18XXh
  oXYH<- P2*P2XYH + P3*P3XYH + P4*P4XYH + P5*P5XYH + P6*P6XYH + P7*P7XYH + P8*P8XYH + P9*P9XYH + P10*P10XYH + P11*P11XYH + P12*P12XYH + P13*P13XYH + P14*P14XYH + P15*P15XYH + P17*P17XYH + P18*P18XYH
  oXYh<- P2*P2XYh + P3*P3XYh + P4*P4XYh + P5*P5XYh + P6*P6XYh + P7*P7XYh + P8*P8XYh + P10*P10XYh + P11*P11XYh + P12*P12XYh + P13*P13XYh + P14*P14XYh + P17*P17XYh + P18*P18XYh
  oXOH<- P2*P2XOH + P3*P3XOH + P4*P4XOH + P5*P5XOH + P6*P6XOH + P7*P7XOH + P8*P8XOH + P9*P9XOH + P10*P10XOH + P11*P11XOH + P12*P12XOH + P13*P13XOH + P14*P14XOH + P15*P15XOH + P17*P17XOH + P18*P18XOH
  oXOh<- P2*P2XOh + P3*P3XOh + P4*P4XOh + P5*P5XOh + P6*P6XOh + P7*P7XOh + P8*P8XOh + P10*P10XOh + P11*P11XOh + P12*P12XOh + P13*P13XOh + P14*P14XOh + P17*P17XOh + P18*P18XOh
  oXOw<- P10*P10XOw + P11*P11XOw + P16*P16XOw
  #all OY and failed conversions are nonviable
  onv <- P2*P2nv + P3*P3nv + P4*P4nv + P5*P5nv + P6*P6nv + P7*P7nv + P8*P8nv + P10*P10nv + P11*P11nv + P12*P12nv + P13*P13nv + P14*P14nv + P17*P17nv + P18*P18nv
  
  #E. Establish matrix of birth probabilities
  pBirth<-matrix(c(oXXw,oXYw,oXXH,oXXh,oXYH,oXYh,oXOH,oXOh,oXOw,onv), nrow=N*N, ncol=9, byrow=F)
  pBirth[is.na(pBirth)] = 0
  
  #F. Adjust offspring probabilities due to fitness costs. GD inds have reduced likelihood of birth
  oXXH<- oXXH*(1-s)
  oXYH<- oXYH*(1-s)
  oXOH<- oXOH*(1-s)
  oXXh<- oXXh*(1-sh)
  oXYh<- oXYh*(1-sh)
  oXOh<- oXOh*(1-sh)
  
  #G. Establish final integer number of offspring of each genotype
  Births <- randBirthYS(Birth = Birthnull, pBirth = pBirth, B = B)
  B1 <- Births[,1]
  B2 <- Births[,2]
  B3 <- Births[,3]
  B4 <- Births[,4]
  B5 <- Births[,5]
  B6 <- Births[,6]
  B7 <- Births[,7]
  B8 <- Births[,8]
  B9 <- Births[,9]
  
  ################################################################################
  #Step 3: Adult dispersal
  
  #Determine number of dispersers
  M1 <- (rbinom(n = mXXw, prob = dA, size = mXXw))
  M2 <- (rbinom(n = mXYw, prob = dA, size = mXYw))
  M3 <- (rbinom(n = mXXH, prob = dA, size = mXXH))
  M4 <- (rbinom(n = mXXh, prob = dA, size = mXXh))
  M5 <- (rbinom(n = mXYH, prob = dA, size = mXYH))
  M6 <- (rbinom(n = mXYh, prob = dA, size = mXYh))
  M7 <- (rbinom(n = mXOH, prob = dA, size = mXOH))
  M8 <- (rbinom(n = mXOh, prob = dA, size = mXOh))
  M9 <- (rbinom(n = mXOw, prob = dA, size = mXOw))
  
  #Returns vector of new immigrants in each patch
  tM1 <- (natMigration2(nMig=M1, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM2 <- (natMigration2(nMig=M2, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM3 <- (natMigration2(nMig=M3, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM4 <- (natMigration2(nMig=M4, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM5<-  (natMigration2(nMig=M5, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM6 <- (natMigration2(nMig=M6, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))   
  tM7 <- (natMigration2(nMig=M7, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM8 <- (natMigration2(nMig=M8, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM9 <- (natMigration2(nMig=M9, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  
  #Update adult pops
  nXXw <- mXXw + tM1 - M1
  nXYw <- mXYw + tM2 - M2
  nXXH <- mXXH + tM3 - M3
  nXXh <- mXXh + tM4 - M4
  nXYH <- mXYH + tM5 - M5
  nXYh <- mXYh + tM6 - M6
  nXOH <- mXOH + tM7 - M7
  nXOh <- mXOh + tM8 - M8
  nXOw <- mXOw + tM9 - M9
  #Ensure population does not go below 1. At least one needed to mate
  nXXw <- ifelse(nXXw<1, 0, nXXw) 
  nXYw <- ifelse(nXYw<1, 0, nXYw) 
  nXXH <- ifelse(nXXH<1, 0, nXXH) 
  nXXh <- ifelse(nXXh<1, 0, nXXh) 
  nXYH <- ifelse(nXYH<1, 0, nXYH) 
  nXYh <- ifelse(nXYh<1, 0, nXYh) 
  nXOH <- ifelse(nXOH<1, 0, nXOH) 
  nXOh <- ifelse(nXOh<1, 0, nXOh)
  nXOw <- ifelse(nXOw<1, 0, nXOw)
  
  ################################################################################
  #Step 4: Juvenile dispersal
  
  #Calc number of dispersers
  M10 <- (rbinom(n = B1, prob = dJ, size = B1))
  M11 <- (rbinom(n = B2, prob = dJ, size = B2))
  M12 <- (rbinom(n = B3, prob = dJ, size = B3))
  M13 <- (rbinom(n = B4, prob = dJ, size = B4))
  M14 <- (rbinom(n = B5, prob = dJ, size = B5))
  M15 <- (rbinom(n = B6, prob = dJ, size = B6))
  M16 <- (rbinom(n = B7, prob = dJ, size = B7))
  M17 <- (rbinom(n = B8, prob = dJ, size = B8))
  M18 <- (rbinom(n = B9, prob = dJ, size = B9))
  
  #Return number of new immigrants in each patch
  tM10 <- (natMigration2(nMig=M10, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM11 <- (natMigration2(nMig=M11, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM12 <- (natMigration2(nMig=M12, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM13 <- (natMigration2(nMig=M13, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM14 <- (natMigration2(nMig=M14, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM15 <- (natMigration2(nMig=M15, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM16 <- (natMigration2(nMig=M16, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM17 <- (natMigration2(nMig=M17, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM18 <- (natMigration2(nMig=M18, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  
  #Update Juvenile populations prior to mortality/maturation
  njXXw <- jXXw + B1 - M10 + tM10
  njXYw <- jXYw + B2 - M11 + tM11
  njXXH <- jXXH + B3 - M12 + tM12
  njXXh <- jXXh + B4 - M13 + tM13
  njXYH <- jXYH + B5 - M14 + tM14
  njXYh <- jXYh + B6 - M15 + tM15
  njXOH <- jXOH + B7 - M16 + tM16
  njXOh <- jXOh + B8 - M17 + tM17
  njXOw <- jXOw + B9 - M18 + tM18
  
  #Setting negative values and NAs to zero
  njXXw <- ifelse(njXXw<0, 0,njXXw)
  njXYw <- ifelse(njXYw<0, 0,njXYw)
  njXXH <- ifelse(njXXH<0, 0,njXXH)
  njXXh <- ifelse(njXXh<0, 0,njXXh)
  njXYH <- ifelse(njXYH<0, 0,njXYH)
  njXYh <- ifelse(njXYh<0, 0,njXYh)
  njXOH <- ifelse(njXOH<0, 0,njXOH)
  njXOh <- ifelse(njXOh<0, 0,njXOh)
  njXOw <- ifelse(njXOw<0, 0,njXOw)
  njXXw[is.na(njXXw)] = 0
  njXYw[is.na(njXYw)] = 0
  njXXH[is.na(njXXH)] = 0
  njXXh[is.na(njXXh)] = 0
  njXYH[is.na(njXYH)] = 0
  njXYh[is.na(njXYh)] = 0
  njXOH[is.na(njXOH)] = 0
  njXOh[is.na(njXOh)] = 0
  njXOw[is.na(njXOw)] = 0
  
  ################################################################################
  #Step 5: Record important results
  
  #Sum births
  totBirth<-sum(B)
  
  #Step 5: Calculate end of time-step totals
  subTotal <- nXXw+nXYw+nXOw+nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXw+njXYw+njXOw+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh
  PopTot <- sum(subTotal) #Global population
  
  GDMO <- nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh
  GDMOTot <- sum(GDMO)   #Global GDMO
  
  #Binary vector indicating if GDMOs are present (1) or not (0) in each patch
  Present <- ifelse(GDMO<1,0,1)
  Present <- as.vector(Present)
  
  #Freq of gd allele in male pop
  GDPatch<- (nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh)/subTotal
  
  #Sex Ratio
  sexRatio<- (nXYw+nXYH+nXYh+njXYw+njXYH+njXYh)/subTotal
  
  
  #Out1: global frequency of gdmo individuals in context of total population
  #      vector of size NL*NL
  GDFrq <- GDMOTot/PopTot #Global Frequency
  GDFrq[is.na(GDFrq)] = 0
  
  #Out2: resistant individual freq
  ResFrq<-0

  #Out3: frequency of cells with at least 1 gdmo adult present (invaded)
  InvFrq <- sum(Present)/(N*N)
  
  #Out4: frequency of patches that have 0 individuals wt or gdmo (extinct)
  Extinct <- ifelse(subTotal<1, 1, 0)
  ExtFrq <- sum(Extinct)/(N*N)
  
  #True Female Mating Rate
  TrueFMR<-FnMated/(mXXw + mXXH + mXXh + mXOH + mXOh + mXOw)
  
  
  #Tracking the time step...maybe not necessary
  ntime.step <-(time.step+1)
  
  dx <- c(nXXw,nXYw,nXOw,nXXH,nXXh,nXYH,nXYh,nXOH,nXOh,njXXw,njXYw,njXOw,njXXH,njXXh,njXYH,njXYh,njXOH,njXOh,
          subTotal,sexRatio, B, GDPatch,
          totBirth,GDMOTot,GDFrq,ResFrq,PopTot,
          InvFrq,ExtFrq,ntime.step) #vector of populations. 
  
  ##combine results into a single vector dx
  list(dx)    
}

#X-Shredder gene drive dynamics. Juvenile dispersal frequency mirrors adult dispersal frequency parameter
Xshred.v2 <- function(t, y, parms, signal){
  
  ################################################################################  
  #Setup: Load Parameters and initital abundances
  N<-NL
  g=parms[[1]]            #Development rate
  u=parms[[2]]            #Mortality rate
  k=parms[[3]]            #Litter size
  v=parms[[4]]            #Female mating rate
  CPop=parms[[5]]         #Carrying capacity of each patch
  Andry=parms[[6]]        #Max male mates per female
  Gyny=parms[[7]]         #Max female mates per male
  kerShape=parms[[8]]     #Dispersal kernel shape
  kerRate=parms[[9]]      #dispersal kernel rate
  dJ=parms[[10]]          #Proportion of juveniles dispersing
  dA=parms[[10]]          #Proportion of adults dispersing
  gdInt=parms[[12]]       #number of gene drive individuals introduced
  gdSupp=parms[[13]]      #release schedule (every 1-10yrs), or once
  intLoc=parms[[14]]      #release all in middle (1) or each ind randomly (2)
  e=parms[[15]]           #prob of HEG1 cutting at homologous locus
  spD=parms[[16]]         #sperm competition disadvantage for drive males (XD)
  Dyr=parms[[18]]         #how many iterations between dispersal events
  
  #Setup: Create vector of N from vector loaded into model
  XY <- y[(0*(N^2) + 1):(1*(N^2))] 
  XD <- y[(1*(N^2) + 1):(2*(N^2))]
  XX <- y[(2*(N^2) + 1):(3*(N^2))]
  jXY <- y[(3*(N^2) + 1):(4*(N^2))]
  jXD <- y[(4*(N^2) + 1):(5*(N^2))]
  jXX <- y[(5*(N^2) + 1):(6*(N^2))]
  time.step<-y[10*(N^2)+10]
  
  ################################################################################  
  #GD Release at specific time steps
  if(gdSupp==11){
    if(time.step==7){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- gdInt #Add specific n of GD ind to central patch 
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
  }
  
  if(gdSupp<11) {
    if(time.step==7){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- gdInt #Always add at the beginning
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
    
    if(time.step>7 & (((time.step-7)/6)%%gdSupp)==0){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- (XD[(N^2/2)+0.5]+gdInt) #Add new GD individuals to central patch
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
  }
  
  ################################################################################
  #Step 1: Mortality and Maturation
  
  #Adult mortality
  D1 <- (rbinom(n = XY, size=XY, prob = u))
  D2 <- (rbinom(n = XD, size=XD, prob = u))
  D3 <- (rbinom(n = XX, size=XX, prob = u))
  
  #Juvenile mortality
  D4 <- (rbinom(n = jXY, size=jXY, prob = u))
  D5 <- (rbinom(n = jXD, size=jXD, prob = u))
  D6 <- (rbinom(n = jXX, size=jXX, prob = u))
  
  #Update juvenile population after mortality
  jXY <- jXY - D4
  jXD <- jXD - D5
  jXX <- jXX - D6
  jXY <- ifelse(jXY<0,0,jXY) 
  jXD <- ifelse(jXD<0,0,jXD) 
  jXX <- ifelse(jXX<0,0,jXX) 
  jXY[is.na(jXY)] = 0
  jXD[is.na(jXD)] = 0
  jXX[is.na(jXX)] = 0
  
  #Juvenile maturation to adult
  J1 <- (rbinom(n = jXY, size=jXY, prob = g))
  J2 <- (rbinom(n = jXD, size=jXD, prob = g))
  J3 <- (rbinom(n = jXX, size=jXX, prob = g))
  
  #Update adult subpopulations
  mXY <- XY + J1 - D1
  mXD <- XD + J2 - D2
  mXX <- XX + J3 - D3
  mXY <- ifelse(mXY<0,0,mXY) 
  mXD <- ifelse(mXD<0,0,mXD) 
  mXX <- ifelse(mXX<0,0,mXX) 
  mXY[is.na(mXY)] = 0
  mXD[is.na(mXD)] = 0
  mXX[is.na(mXX)] = 0
  
  #Update juvenile population after maturation
  jXY <- jXY - J1
  jXD <- jXD - J2
  jXX <- jXX - J3
  jXY <- ifelse(jXY<0,0,jXY) 
  jXD <- ifelse(jXD<0,0,jXD) 
  jXX <- ifelse(jXX<0,0,jXX) 
  jXY[is.na(jXY)] = 0
  jXD[is.na(jXD)] = 0
  jXX[is.na(jXX)] = 0
  
  ################################################################################
  #Step 2: Adult mating
  
  #A. Mating population established (mating males "mM" and mating females "mF")
  mM <- mXY + mXD   #Male mating population
  mF <- v*(mXX) #Female Mating Population 
  
  #Number mated females given max number of mates per males and females
  #Requires polyandry probabilities (likelihood of n mates for females for 1:max)
  
  if(Andry==1){
    FnMated<-mapply(minCheck, round(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq)), mF)
    FnMated[is.na(FnMated)] = 0 
  }else{
    FnMated<-mapply(minCheck, round(colSums(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq))), mF)
    FnMated[is.na(FnMated)] = 0 
  }
  
  #B. Total number of births in each patch w density dependence, given avg litter size "k"
  sumLit<-(randpoisLitter(FnMated, k, FnMated))  #assumes min=2, max=9
  B<- sumLit*(1-((mXY + mXD + mXX)/CPop))
  B<- ifelse(B<0, 0, B)
  B <- (randpois(dsubpop = B, lambda = 1, subpop = B))
  
  #C. Mating fractions among all adult genotype pairs, assuming random mating within spatial patches
  allmate<- (mXX * (mXY + mXD))
  P1<- (mXX * mXY) /allmate
  P2<- (mXX * mXD) /allmate
  
  #D. calculate probability of offspring genotypes
  oXY<- P1*P1XY
  oXD<- P2*P2XD
  oXX<- P1*P1XX + P2*P2XX
  
  #E. Adjust offspring probabilities due to fitness costs. GD inds have reduced likelihood of birth
  if(length(PolyandrySeq)>1){
    oXD <- oXD-(oXD*(1-spD)*(sum(PolyandryProbs[-1])))
  }
  
  #F. Establish matrix of mating probabilities for randBirth2 function
  pBirth <-matrix(c(oXY,oXD,oXX),nrow=N*N,ncol=3,byrow=FALSE) #Male Bias. 
  
  #G. Establish final integer number of offpsring of each genotype
  Births <- randBirthXS(Birth = Birthnull, pBirth = pBirth, B = B)
  B1 <- Births[,1]
  B2 <- Births[,2]
  B3 <- Births[,3]
  
  
  ################################################################################
  #Step 3: Adult dispersal

  if((time.step-8)%%Dyr==0){
    
    M1 <- (rbinom(n = mXY, size=mXY, prob = dA))
    M2 <- (rbinom(n = mXD, size=mXD, prob = dA))
    M3 <- (rbinom(n = mXX, size=mXX, prob = dA))
    
    #Returns vector of new immigrants in each patch
    #Needs shape and rate parameters
    #Needs pSize parameter for informing range of potential destination patches
    tM1 <- (natMigration2(nMig=M1, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM2 <- (natMigration2(nMig=M2, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM3 <- (natMigration2(nMig=M3, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    
    #Update adult pops
    nXY <- mXY + tM1 - M1
    nXD <- mXD + tM2 - M2
    nXX <- mXX + tM3 - M3
    nXY <- ifelse(nXY<1, 0, nXY) #Ensure population does not go below 1. At least one needed to mate
    nXD <- ifelse(nXD<1, 0, nXD) 
    nXX <- ifelse(nXX<1, 0, nXX) 
  } else{
    nXY <- mXY
    nXD <- mXD
    nXX <- mXX
    nXY <- ifelse(nXY<1, 0, nXY) #Ensure population does not go below 1. At least one needed to mate
    nXD <- ifelse(nXD<1, 0, nXD) 
    nXX <- ifelse(nXX<1, 0, nXX) 
    
  }
  
  ##################################################################################
  #Step 4: Juvenile Dispersal
  
  if((time.step-8)%%Dyr==0){
    
    #Calc number of dispersers
    M4 <- (rbinom(n = B1, size=B1, prob = dJ))
    M5 <- (rbinom(n = B2, size=B2, prob = dJ))
    M6 <- (rbinom(n = B3, size=B3, prob = dJ))
    
    #Return number of new immigrants in each patch
    tM4 <- (natMigration2(nMig=M4, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM5 <- (natMigration2(nMig=M5, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM6 <- (natMigration2(nMig=M6, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    
    #Update Juvenile populations from births and dispersal
    njXY <- jXY + B1 - M4 + tM4 
    njXD <- jXD + B2 - M5 + tM5
    njXX <- jXX + B3 - M6 + tM6
    
    njXY <- ifelse(njXY<0, 0,njXY)
    njXD <- ifelse(njXD<0, 0,njXD)
    njXX <- ifelse(njXX<0, 0,njXX)
    njXY[is.na(njXY)] = 0
    njXD[is.na(njXD)] = 0
    njXX[is.na(njXX)] = 0
  } else {
    njXY <- jXY + B1 
    njXD <- jXD + B2
    njXX <- jXX + B3
    
    njXY <- ifelse(njXY<0, 0,njXY)
    njXD <- ifelse(njXD<0, 0,njXD)
    njXX <- ifelse(njXX<0, 0,njXX)
    njXY[is.na(njXY)] = 0
    njXD[is.na(njXD)] = 0
    njXX[is.na(njXX)] = 0
  }
  ################################################################################
  #Step 5: Record important results
  
  #####
  #Sum births
  totBirth<-sum(B)
  
  #Step 5: Calculate end of time-step totals
  patchTotal <- nXY+nXD+nXX+njXY+njXD+njXX
  PopTot <- sum(patchTotal) #Global population
  
  GDMO <- nXD  #not including nRD because dead end individuals that wont drive themselves
  GDMOTot <- sum(GDMO)   #Global GDMO
  
  #Binary vector indicating if GDMOs are present (1) or not (0) in each patch
  Present <- ifelse(GDMO<1, 0,1)
  Present <- as.vector(Present)
  
  #Freq of gd allele in male pop
  GDMales<- (nXD)/(nXY+nXD)
  
  #genotype total for animation
  fqnXD<-sum(nXD+njXD)/PopTot
  fqnXY<-sum(nXY+njXY)/PopTot
  fqnXX<-sum(nXX+njXX)/PopTot
  
  #Sex Ratio in each patch (>0.5 is skewed male, <0.5 is skewed female)
  sexRatio<- (nXY+nXD+njXY+njXD) / patchTotal
  
  #Number of births in each patch
  #already contained in vector B
  
  #Out1: global frequency of gdmo individuals in context of total population
  #      vector of size NL*NL
  GDFrq <- GDMOTot/PopTot #Global Frequency
  GDFrq[is.na(GDFrq)] = 0
  
  #Out3: frequency of cells with at least 1 gdmo adult present (invaded)
  InvFrq <- sum(Present)/(N*N)
  
  #Out4: frequency of patches that have 0 individuals wt or gdmo (extinct)
  Extinct <- ifelse(patchTotal<1, 1, 0)
  ExtFrq <- sum(Extinct)/(N*N)
  
  #Tracking the time step...maybe not necessary
  ntime.step <-(time.step+1)
  
  dx <- c(nXY,nXD,nXX,njXY,njXD,njXX,
          patchTotal, sexRatio, B, GDMales,
          totBirth, GDMOTot, GDFrq, PopTot,
          InvFrq, ExtFrq, fqnXD, fqnXY, fqnXX, ntime.step) #vector of populations.
  
  #######################################
  ##combine results into a single vector dx
  ######################################
  list(dx)  
}

#Y-Shredder gene drive dynamics. Juvenile dispersal frequency mirrors adult dispersal frequency parameter
Yshred.v2 <- function(t, y, parms, signal){
  
  ################################################################################  
  #Setup: Load Parameters and initital abundances
  
  N<-NL
  g=parms[[1]]            #Development rate
  u=parms[[2]]            #Mortality rate
  k=parms[[3]]            #Litter size
  v=parms[[4]]            #Female mating rate
  CPop=parms[[5]]         #Carrying capacity of each patch
  Andry=parms[[6]]        #Max male mates per female
  Gyny=parms[[7]]         #Max female mates per male
  kerShape=parms[[8]]     #Dispersal kernel shape
  kerRate=parms[[9]]      #dispersal kernel rate
  dJ=parms[[10]]          #Proportion of juveniles dispersing
  dA=parms[[10]]          #Proportion of adults dispersing
  gdInt=parms[[12]]       #number of gene drive individuals introduced
  gdSupp=parms[[13]]      #release schedule (every 1-10yrs), or once
  intLoc=parms[[14]]      #release all in middle (1) or each ind randomly (2)
  Pc=parms[[15]]          #prob of HEG1 cutting at homologous locus
  Pn=parms[[16]]          #prob of NHEJ after HEG1 cut (i.e. successful conversion)
  Py=parms[[17]]          #prob of HEG2 cutting at Y chromosome
  Xpass=parms[[18]]       #prob of XO female passing X chrom
  s=parms[[19]]           #Fitness cost for homozygous GDMO
  sh=parms[[20]]          #Fitness cost for heterozygous GDMO
  
  #Setup: Create vector of N from vector loaded into model
  XXw <- y[(0*(N^2) + 1): (1*(N^2))] 
  XYw <- y[(1*(N^2) + 1): (2*(N^2))]
  XOw <- y[(2*(N^2) + 1): (3*(N^2))]
  XXH <- y[(3*(N^2) + 1): (4*(N^2))]
  XXh <- y[(4*(N^2) + 1): (5*(N^2))]
  XYH <- y[(5*(N^2) + 1): (6*(N^2))]
  XYh <- y[(6*(N^2) + 1): (7*(N^2))]
  XOH <- y[(7*(N^2) + 1): (8*(N^2))]
  XOh <- y[(8*(N^2) + 1): (9*(N^2))]
  jXXw <- y[(9*(N^2) + 1): (10*(N^2))]
  jXYw <- y[(10*(N^2) + 1): (11*(N^2))]
  jXOw <- y[(11*(N^2) + 1): (12*(N^2))]
  jXXH <- y[(12*(N^2) + 1): (13*(N^2))]
  jXXh <- y[(13*(N^2) + 1): (14*(N^2))]
  jXYH <- y[(14*(N^2) + 1): (15*(N^2))]
  jXYh <- y[(15*(N^2) + 1): (16*(N^2))]
  jXOH <- y[(16*(N^2) + 1): (17*(N^2))]
  jXOh <- y[(17*(N^2) + 1): (18*(N^2))]
  time.step<-y[((22*(N^2))+8)]
  
  ################################################################################  
  #GD Release at specific time steps
  
  #Scheduled release of GD individuals (heterozygous females), once here, successively below
  if(gdSupp==11){
    if(time.step==7){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- gdInt #Add specific n of GD ind to central patch 
      }
      if(intLoc==2){              #Random placement
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
  }
  #First release when supplemental releases scheduled
  if(gdSupp<11) {
    if(time.step==7){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- gdInt #Always add at the beginning
      }
      if(intLoc==2){              #Random placement
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
    #Supplemental release every X years
    if(time.step>7 & (((time.step-7)/6)%%gdSupp)==0){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- (XXH[(N^2/2)+0.5]+gdInt) #Add new GD individuals to central patch
      }
      if(intLoc==2){
        for(i in 1:gdInt){              #Random placement
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
  }
  
  ################################################################################
  #Step 1: Mortality and Maturation
  
  #Determine number of adult mortalities
  D1 <- (rbinom(n = XXw, prob = u, size = XXw))
  D2 <- (rbinom(n = XYw, prob = u, size = XYw))
  D3 <- (rbinom(n = XXH, prob = u, size = XXH))
  D4 <- (rbinom(n = XXh, prob = u, size = XXh))
  D5 <- (rbinom(n = XYH, prob = u, size = XYH))
  D6 <- (rbinom(n = XYh, prob = u, size = XYh))
  D7 <- (rbinom(n = XOH, prob = u, size = XOH))
  D8 <- (rbinom(n = XOh, prob = u, size = XOh))
  D9 <- (rbinom(n = XOw, prob = u, size = XOw))
  
  #Determine number of juvenile mortalities
  D10 <- (rbinom(n = jXXw, prob = u, size = jXXw))
  D11 <- (rbinom(n = jXYw, prob = u, size = jXYw))
  D12 <- (rbinom(n = jXXH, prob = u, size = jXXH))
  D13 <- (rbinom(n = jXXh, prob = u, size = jXXh))
  D14 <- (rbinom(n = jXYH, prob = u, size = jXYH))
  D15 <- (rbinom(n = jXYh, prob = u, size = jXYh))
  D16 <- (rbinom(n = jXOH, prob = u, size = jXOH))
  D17 <- (rbinom(n = jXOh, prob = u, size = jXOh))
  D18 <- (rbinom(n = jXOw, prob = u, size = jXOw))
  
  #Update juvenile population after mortality
  jXXw <- jXXw - D10
  jXYw <- jXYw - D11
  jXXH <- jXXH - D12
  jXXh <- jXXh - D13
  jXYH <- jXYH - D14
  jXYh <- jXYh - D15
  jXOH <- jXOH - D16
  jXOh <- jXOh - D17
  jXOw <- jXOw - D18
  jXXw <- ifelse(jXXw<0,0,jXXw) 
  jXYw <- ifelse(jXYw<0,0,jXYw) 
  jXXH <- ifelse(jXXH<0,0,jXXH) 
  jXXh <- ifelse(jXXh<0,0,jXXh) 
  jXYH <- ifelse(jXYH<0,0,jXYH) 
  jXYh <- ifelse(jXYh<0,0,jXYh) 
  jXOH <- ifelse(jXOH<0,0,jXOH) 
  jXOh <- ifelse(jXOh<0,0,jXOh) 
  jXOw <- ifelse(jXOw<0,0,jXOw) 
  jXXw[is.na(jXXw)] = 0
  jXYw[is.na(jXYw)] = 0
  jXXH[is.na(jXXH)] = 0
  jXXh[is.na(jXXh)] = 0
  jXYH[is.na(jXYH)] = 0
  jXYh[is.na(jXYh)] = 0
  jXOH[is.na(jXOH)] = 0
  jXOh[is.na(jXOh)] = 0
  jXOw[is.na(jXOw)] = 0
  
  #Poisson of maturation
  J1 <- (rbinom(n = jXXw, prob = g, size = jXXw))
  J2 <- (rbinom(n = jXYw, prob = g, size = jXYw))
  J3 <- (rbinom(n = jXXH, prob = g, size = jXXH))
  J4 <- (rbinom(n = jXXh, prob = g, size = jXXh))
  J5 <- (rbinom(n = jXYH, prob = g, size = jXYH))
  J6 <- (rbinom(n = jXYh, prob = g, size = jXYh))
  J7 <- (rbinom(n = jXOH, prob = g, size = jXOH))
  J8 <- (rbinom(n = jXOh, prob = g, size = jXOh))
  J9 <- (rbinom(n = jXOw, prob = g, size = jXOw))
  
  #Update adult subpopulations
  mXXw <- XXw + J1 - D1
  mXYw <- XYw + J2 - D2
  mXXH <- XXH + J3 - D3
  mXXh <- XXh + J4 - D4
  mXYH <- XYH + J5 - D5
  mXYh <- XYh + J6 - D6
  mXOH <- XOH + J7 - D7
  mXOh <- XOh + J8 - D8
  mXOw <- XOw + J9 - D9
  mXXw <- ifelse(mXXw<0,0,mXXw) 
  mXYw <- ifelse(mXYw<0,0,mXYw) 
  mXXH <- ifelse(mXXH<0,0,mXXH) 
  mXXh <- ifelse(mXXh<0,0,mXXh) 
  mXYH <- ifelse(mXYH<0,0,mXYH) 
  mXYh <- ifelse(mXYh<0,0,mXYh)
  mXOH <- ifelse(mXOH<0,0,mXOH)
  mXOh <- ifelse(mXOh<0,0,mXOh)
  mXOw <- ifelse(mXOw<0,0,mXOw)
  mXXw[is.na(mXXw)] = 0
  mXYw[is.na(mXYw)] = 0
  mXXH[is.na(mXXH)] = 0
  mXXh[is.na(mXXh)] = 0
  mXYH[is.na(mXYH)] = 0
  mXYh[is.na(mXYh)] = 0  
  mXOH[is.na(mXOH)] = 0 
  mXOh[is.na(mXOh)] = 0 
  mXOw[is.na(mXOw)] = 0 
  
  #Update juvenile population after maturation
  jXXw <- jXXw - J1
  jXYw <- jXYw - J2
  jXXH <- jXXH - J3
  jXXh <- jXXh - J4
  jXYH <- jXYH - J5
  jXYh <- jXYh - J6
  jXOH <- jXOH - J7
  jXOh <- jXOh - J8
  jXOw <- jXOw - J9
  jXXw <- ifelse(jXXw<0,0,jXXw) 
  jXYw <- ifelse(jXYw<0,0,jXYw) 
  jXXH <- ifelse(jXXH<0,0,jXXH) 
  jXXh <- ifelse(jXXh<0,0,jXXh) 
  jXYH <- ifelse(jXYH<0,0,jXYH) 
  jXYh <- ifelse(jXYh<0,0,jXYh) 
  jXOH <- ifelse(jXOH<0,0,jXOH) 
  jXOh <- ifelse(jXOh<0,0,jXOh) 
  jXOw <- ifelse(jXOw<0,0,jXOw) 
  jXXw[is.na(jXXw)] = 0
  jXYw[is.na(jXYw)] = 0
  jXXH[is.na(jXXH)] = 0
  jXXh[is.na(jXXh)] = 0
  jXYH[is.na(jXYH)] = 0
  jXYh[is.na(jXYh)] = 0
  jXOH[is.na(jXOH)] = 0
  jXOh[is.na(jXOh)] = 0
  jXOw[is.na(jXOw)] = 0
  
  ################################################################################
  #Step 2: Adult mating
  
  #A. Mating population established (mating males "mM" and mating females "mF")
  mM <- mXYw + mXYH + mXYh                    #Male mating population
  mF <- v*(mXXw + mXXH + mXXh + mXOH + mXOh + mXOw)  #Female Mating Population 
  
  # Determine number mated females given max number of mates per males and females
  #Requires polyandry probabilities (likelihood of n mates for females for 1:max)
  if(Andry==1){
    FnMated<-mapply(minCheck, round(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq)), mF)
    FnMated[is.na(FnMated)] = 0 
  }else{
    FnMated<-mapply(minCheck, round(colSums(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq))), mF)
    FnMated[is.na(FnMated)] = 0 
  }
  
  #B. Total number of births in each patch w density dependence, given avg litter size "k"
  sumLit<-(randpoisLitter(FnMated, k, FnMated))  #assumes min=2, max=9
  B<- sumLit*(1-((mXXw + mXYw + mXXH + mXXh + mXYH + mXYh + mXOH + mXOh + mXOw)/CPop))
  B<- ifelse(B<0, 0,B)
  B <- (randpois(dsubpop = B, lambda = 1, subpop = B))
  
  #C. Mating fractions among all adult genotype pairs, assuming random mating within spatial patches
  allmate<- (mXXw * (mXYw + mXYH +mXYh)) + (mXXH * (mXYw + mXYH +mXYh)) + (mXXh * (mXYw + mXYH +mXYh)) + (mXOH * (mXYw + mXYH +mXYh)) + (mXOh * (mXYw + mXYH +mXYh)) + (mXOw * (mXYw + mXYH +mXYh))
  P1<-(mXXw * mXYw) / allmate
  P2<-(mXXw * mXYh) / allmate
  P3<-(mXXw * mXYH) / allmate
  P4<-(mXXh * mXYw) / allmate
  P5<-(mXXh * mXYh) / allmate
  P6<-(mXXh * mXYH) / allmate
  P7<-(mXXH * mXYw) / allmate
  P8<-(mXXH * mXYh) / allmate
  P9<-(mXXH * mXYH) / allmate
  P10<-(mXOh * mXYw) / allmate
  P11<-(mXOh * mXYh) / allmate
  P12<-(mXOh * mXYH) / allmate
  P13<-(mXOH * mXYw) / allmate
  P14<-(mXOH * mXYh) / allmate
  P15<-(mXOH * mXYH) / allmate
  P16<-(mXOw * mXYw) / allmate
  P17<-(mXOw * mXYH) / allmate
  P18<-(mXOw * mXYh) / allmate
  
  #D. Offspring probs sum of mating fraction and inheritance calculation for all possible parent combos that could produce offspring geno
  oXXw<- P1*P1XXw + P2*P2XXw + P4*P4XXw + P5*P5XXw + P10*P10XXw + P11*P11XXw + P15*P15XXw + P18*P18XXw
  oXYw<- P1*P1XYw + P2*P2XYw + P4*P4XYw + P5*P5XYw + P10*P10XYw + P11*P11XYw + P15*P15XYw + P18*P18XYw
  oXXH<- P2*P2XXH + P3*P3XXH + P4*P4XXH + P5*P5XXH + P6*P6XXH + P7*P7XXH + P8*P8XXH + P9*P9XXH + P10*P10XXH + P11*P11XXH + P12*P12XXH + P13*P13XXH + P14*P14XXH + P15*P15XXH + P17*P17XXH + P18*P18XXH
  oXXh<- P2*P2XXh + P3*P3XXh + P4*P4XXh + P5*P5XXh + P6*P6XXh + P7*P7XXh + P8*P8XXh + P10*P10XXh + P11*P11XXh + P12*P12XXh + P13*P13XXh + P14*P14XXh + P17*P17XXh + P18*P18XXh
  oXYH<- P2*P2XYH + P3*P3XYH + P4*P4XYH + P5*P5XYH + P6*P6XYH + P7*P7XYH + P8*P8XYH + P9*P9XYH + P10*P10XYH + P11*P11XYH + P12*P12XYH + P13*P13XYH + P14*P14XYH + P15*P15XYH + P17*P17XYH + P18*P18XYH
  oXYh<- P2*P2XYh + P3*P3XYh + P4*P4XYh + P5*P5XYh + P6*P6XYh + P7*P7XYh + P8*P8XYh + P10*P10XYh + P11*P11XYh + P12*P12XYh + P13*P13XYh + P14*P14XYh + P17*P17XYh + P18*P18XYh
  oXOH<- P2*P2XOH + P3*P3XOH + P4*P4XOH + P5*P5XOH + P6*P6XOH + P7*P7XOH + P8*P8XOH + P9*P9XOH + P10*P10XOH + P11*P11XOH + P12*P12XOH + P13*P13XOH + P14*P14XOH + P15*P15XOH + P17*P17XOH + P18*P18XOH
  oXOh<- P2*P2XOh + P3*P3XOh + P4*P4XOh + P5*P5XOh + P6*P6XOh + P7*P7XOh + P8*P8XOh + P10*P10XOh + P11*P11XOh + P12*P12XOh + P13*P13XOh + P14*P14XOh + P17*P17XOh + P18*P18XOh
  oXOw<- P10*P10XOw + P11*P11XOw + P16*P16XOw
  #all OY and failed conversions are nonviable
  onv <- P2*P2nv + P3*P3nv + P4*P4nv + P5*P5nv + P6*P6nv + P7*P7nv + P8*P8nv + P10*P10nv + P11*P11nv + P12*P12nv + P13*P13nv + P14*P14nv + P17*P17nv + P18*P18nv
  
  #E. Establish matrix of birth probabilities
  pBirth<-matrix(c(oXXw,oXYw,oXXH,oXXh,oXYH,oXYh,oXOH,oXOh,oXOw,onv), nrow=N*N, ncol=9, byrow=F)
  pBirth[is.na(pBirth)] = 0
  
  #F. Adjust offspring probabilities due to fitness costs. GD inds have reduced likelihood of birth
  oXXH<- oXXH*(1-s)
  oXYH<- oXYH*(1-s)
  oXOH<- oXOH*(1-s)
  oXXh<- oXXh*(1-sh)
  oXYh<- oXYh*(1-sh)
  oXOh<- oXOh*(1-sh)
  
  #G. Establish final integer number of offspring of each genotype
  Births <- randBirthYS(Birth = Birthnull, pBirth = pBirth, B = B)
  B1 <- Births[,1]
  B2 <- Births[,2]
  B3 <- Births[,3]
  B4 <- Births[,4]
  B5 <- Births[,5]
  B6 <- Births[,6]
  B7 <- Births[,7]
  B8 <- Births[,8]
  B9 <- Births[,9]
  
  ################################################################################
  #Step 3: Adult dispersal
  
  #Determine number of dispersers
  M1 <- (rbinom(n = mXXw, prob = dA, size = mXXw))
  M2 <- (rbinom(n = mXYw, prob = dA, size = mXYw))
  M3 <- (rbinom(n = mXXH, prob = dA, size = mXXH))
  M4 <- (rbinom(n = mXXh, prob = dA, size = mXXh))
  M5 <- (rbinom(n = mXYH, prob = dA, size = mXYH))
  M6 <- (rbinom(n = mXYh, prob = dA, size = mXYh))
  M7 <- (rbinom(n = mXOH, prob = dA, size = mXOH))
  M8 <- (rbinom(n = mXOh, prob = dA, size = mXOh))
  M9 <- (rbinom(n = mXOw, prob = dA, size = mXOw))
  
  #Returns vector of new immigrants in each patch
  tM1 <- (natMigration2(nMig=M1, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM2 <- (natMigration2(nMig=M2, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM3 <- (natMigration2(nMig=M3, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM4 <- (natMigration2(nMig=M4, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM5<-  (natMigration2(nMig=M5, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM6 <- (natMigration2(nMig=M6, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))   
  tM7 <- (natMigration2(nMig=M7, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM8 <- (natMigration2(nMig=M8, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM9 <- (natMigration2(nMig=M9, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  
  #Update adult pops
  nXXw <- mXXw + tM1 - M1
  nXYw <- mXYw + tM2 - M2
  nXXH <- mXXH + tM3 - M3
  nXXh <- mXXh + tM4 - M4
  nXYH <- mXYH + tM5 - M5
  nXYh <- mXYh + tM6 - M6
  nXOH <- mXOH + tM7 - M7
  nXOh <- mXOh + tM8 - M8
  nXOw <- mXOw + tM9 - M9
  #Ensure population does not go below 1. At least one needed to mate
  nXXw <- ifelse(nXXw<1, 0, nXXw) 
  nXYw <- ifelse(nXYw<1, 0, nXYw) 
  nXXH <- ifelse(nXXH<1, 0, nXXH) 
  nXXh <- ifelse(nXXh<1, 0, nXXh) 
  nXYH <- ifelse(nXYH<1, 0, nXYH) 
  nXYh <- ifelse(nXYh<1, 0, nXYh) 
  nXOH <- ifelse(nXOH<1, 0, nXOH) 
  nXOh <- ifelse(nXOh<1, 0, nXOh)
  nXOw <- ifelse(nXOw<1, 0, nXOw)
  
  ################################################################################
  #Step 4: Juvenile dispersal
  
  #Calc number of dispersers
  M10 <- (rbinom(n = B1, prob = dJ, size = B1))
  M11 <- (rbinom(n = B2, prob = dJ, size = B2))
  M12 <- (rbinom(n = B3, prob = dJ, size = B3))
  M13 <- (rbinom(n = B4, prob = dJ, size = B4))
  M14 <- (rbinom(n = B5, prob = dJ, size = B5))
  M15 <- (rbinom(n = B6, prob = dJ, size = B6))
  M16 <- (rbinom(n = B7, prob = dJ, size = B7))
  M17 <- (rbinom(n = B8, prob = dJ, size = B8))
  M18 <- (rbinom(n = B9, prob = dJ, size = B9))
  
  #Return number of new immigrants in each patch
  tM10 <- (natMigration2(nMig=M10, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM11 <- (natMigration2(nMig=M11, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM12 <- (natMigration2(nMig=M12, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM13 <- (natMigration2(nMig=M13, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM14 <- (natMigration2(nMig=M14, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM15 <- (natMigration2(nMig=M15, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM16 <- (natMigration2(nMig=M16, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM17 <- (natMigration2(nMig=M17, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM18 <- (natMigration2(nMig=M18, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  
  #Update Juvenile populations prior to mortality/maturation
  njXXw <- jXXw + B1 - M10 + tM10
  njXYw <- jXYw + B2 - M11 + tM11
  njXXH <- jXXH + B3 - M12 + tM12
  njXXh <- jXXh + B4 - M13 + tM13
  njXYH <- jXYH + B5 - M14 + tM14
  njXYh <- jXYh + B6 - M15 + tM15
  njXOH <- jXOH + B7 - M16 + tM16
  njXOh <- jXOh + B8 - M17 + tM17
  njXOw <- jXOw + B9 - M18 + tM18
  
  #Setting negative values and NAs to zero
  njXXw <- ifelse(njXXw<0, 0,njXXw)
  njXYw <- ifelse(njXYw<0, 0,njXYw)
  njXXH <- ifelse(njXXH<0, 0,njXXH)
  njXXh <- ifelse(njXXh<0, 0,njXXh)
  njXYH <- ifelse(njXYH<0, 0,njXYH)
  njXYh <- ifelse(njXYh<0, 0,njXYh)
  njXOH <- ifelse(njXOH<0, 0,njXOH)
  njXOh <- ifelse(njXOh<0, 0,njXOh)
  njXOw <- ifelse(njXOw<0, 0,njXOw)
  njXXw[is.na(njXXw)] = 0
  njXYw[is.na(njXYw)] = 0
  njXXH[is.na(njXXH)] = 0
  njXXh[is.na(njXXh)] = 0
  njXYH[is.na(njXYH)] = 0
  njXYh[is.na(njXYh)] = 0
  njXOH[is.na(njXOH)] = 0
  njXOh[is.na(njXOh)] = 0
  njXOw[is.na(njXOw)] = 0
  
  ################################################################################
  #Step 5: Record important results
  
  #Sum births
  totBirth<-sum(B)
  
  #Step 5: Calculate end of time-step totals
  subTotal <- nXXw+nXYw+nXOw+nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXw+njXYw+njXOw+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh
  PopTot <- sum(subTotal) #Global population
  
  GDMO <- nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh
  GDMOTot <- sum(GDMO)   #Global GDMO
  
  #Binary vector indicating if GDMOs are present (1) or not (0) in each patch
  Present <- ifelse(GDMO<1,0,1)
  Present <- as.vector(Present)
  
  #Freq of gd allele in male pop
  GDPatch<- (nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh)/subTotal
  
  #Sex Ratio
  sexRatio<- (nXYw+nXYH+nXYh+njXYw+njXYH+njXYh)/subTotal
  
  
  #Out1: global frequency of gdmo individuals in context of total population
  #      vector of size NL*NL
  GDFrq <- GDMOTot/PopTot #Global Frequency
  GDFrq[is.na(GDFrq)] = 0
  
  #Out2: resistant individual freq
  ResFrq<-0
  
  #Out2: Global abundance over time, single value per time step
  # just use "PopTot"...don't need to make new object
  
  #Out3: frequency of cells with at least 1 gdmo adult present (invaded)
  InvFrq <- sum(Present)/(N*N)
  
  #Out4: frequency of patches that have 0 individuals wt or gdmo (extinct)
  Extinct <- ifelse(subTotal<1, 1, 0)
  ExtFrq <- sum(Extinct)/(N*N)
  
  #True Female Mating Rate
  TrueFMR<-FnMated/(mXXw + mXXH + mXXh + mXOH + mXOh + mXOw)
  
  
  #Tracking the time step...maybe not necessary
  ntime.step <-(time.step+1)
  
  dx <- c(nXXw,nXYw,nXOw,nXXH,nXXh,nXYH,nXYh,nXOH,nXOh,njXXw,njXYw,njXOw,njXXH,njXXh,njXYH,njXYh,njXOH,njXOh,
          subTotal,sexRatio, B, GDPatch,
          totBirth,GDMOTot,GDFrq,ResFrq,PopTot,
          InvFrq,ExtFrq,ntime.step) #vector of populations. 
  
  ##combine results into a single vector dx
  list(dx)    
}

#X-shredder gene drive dynamics. Juvenile dispersal frequency mirrors adult dispersal frequency parameter. Tracks chasing dynamics
Xshred.v2.chase <- function(t, y, parms, signal){
  
  ################################################################################  
  #Setup: Load Parameters and initital abundances
  N<-NL
  g=parms[[1]]            #Development rate
  u=parms[[2]]            #Mortality rate
  k=parms[[3]]            #Litter size
  v=parms[[4]]            #Female mating rate
  CPop=parms[[5]]         #Carrying capacity of each patch
  Andry=parms[[6]]        #Max male mates per female
  Gyny=parms[[7]]         #Max female mates per male
  kerShape=parms[[8]]     #Dispersal kernel shape
  kerRate=parms[[9]]      #dispersal kernel rate
  dJ=parms[[10]]          #Proportion of juveniles dispersing
  dA=parms[[10]]          #Proportion of adults dispersing
  gdInt=parms[[12]]       #number of gene drive individuals introduced
  gdSupp=parms[[13]]      #release schedule (every 1-10yrs), or once
  intLoc=parms[[14]]      #release all in middle (1) or each ind randomly (2)
  e=parms[[15]]           #prob of HEG1 cutting at homologous locus
  spD=parms[[16]]         #sperm competition disadvantage for drive males (XD)
  Dyr=parms[[18]]         #how many iterations between dispersal events
  
  #Setup: Create vector of N from vector loaded into model
  XY <- y[(0*(N^2) + 1):(1*(N^2))] 
  XD <- y[(1*(N^2) + 1):(2*(N^2))]
  XX <- y[(2*(N^2) + 1):(3*(N^2))]
  jXY <- y[(3*(N^2) + 1):(4*(N^2))]
  jXD <- y[(4*(N^2) + 1):(5*(N^2))]
  jXX <- y[(5*(N^2) + 1):(6*(N^2))]
  time.step<-y[10*(N^2)+13]
  
  ################################################################################  
  #GD Release at specific time steps
  if(gdSupp==11){
    if(time.step==7){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- gdInt #Add specific n of GD ind to central patch 
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
  }
  
  if(gdSupp<11) {
    if(time.step==7){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- gdInt #Always add at the beginning
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
    
    if(time.step>7 & (((time.step-7)/6)%%gdSupp)==0){
      if(intLoc==1){
        XD[(N^2/2)+0.5] <- (XD[(N^2/2)+0.5]+gdInt) #Add new GD individuals to central patch
      }
      if(intLoc==2){
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XD[loc]<- (XD[loc]+1)
        }
      }
    }
  }
  
  ################################################################################
  #Step 1: Mortality and Maturation
  
  #Adult mortality
  D1 <- (rbinom(n = XY, size=XY, prob = u))
  D2 <- (rbinom(n = XD, size=XD, prob = u))
  D3 <- (rbinom(n = XX, size=XX, prob = u))
  
  #Juvenile mortality
  D4 <- (rbinom(n = jXY, size=jXY, prob = u))
  D5 <- (rbinom(n = jXD, size=jXD, prob = u))
  D6 <- (rbinom(n = jXX, size=jXX, prob = u))
  
  #Update juvenile population after mortality
  jXY <- jXY - D4
  jXD <- jXD - D5
  jXX <- jXX - D6
  jXY <- ifelse(jXY<0,0,jXY) 
  jXD <- ifelse(jXD<0,0,jXD) 
  jXX <- ifelse(jXX<0,0,jXX) 
  jXY[is.na(jXY)] = 0
  jXD[is.na(jXD)] = 0
  jXX[is.na(jXX)] = 0
  
  #Juvenile maturation to adult
  J1 <- (rbinom(n = jXY, size=jXY, prob = g))
  J2 <- (rbinom(n = jXD, size=jXD, prob = g))
  J3 <- (rbinom(n = jXX, size=jXX, prob = g))
  
  #Update adult subpopulations
  mXY <- XY + J1 - D1
  mXD <- XD + J2 - D2
  mXX <- XX + J3 - D3
  mXY <- ifelse(mXY<0,0,mXY) 
  mXD <- ifelse(mXD<0,0,mXD) 
  mXX <- ifelse(mXX<0,0,mXX) 
  mXY[is.na(mXY)] = 0
  mXD[is.na(mXD)] = 0
  mXX[is.na(mXX)] = 0
  
  #Update juvenile population after maturation
  jXY <- jXY - J1
  jXD <- jXD - J2
  jXX <- jXX - J3
  jXY <- ifelse(jXY<0,0,jXY) 
  jXD <- ifelse(jXD<0,0,jXD) 
  jXX <- ifelse(jXX<0,0,jXX) 
  jXY[is.na(jXY)] = 0
  jXD[is.na(jXD)] = 0
  jXX[is.na(jXX)] = 0
  
  ################################################################################
  #Step 2: Adult mating
  
  #A. Mating population established (mating males "mM" and mating females "mF")
  mM <- mXY + mXD   #Male mating population
  mF <- v*(mXX) #Female Mating Population 
  
  #Number mated females given max number of mates per males and females
  #Requires polyandry probabilities (likelihood of n mates for females for 1:max)
  
  if(Andry==1){
    FnMated<-mapply(minCheck, round(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq)), mF)
    FnMated[is.na(FnMated)] = 0 
  }else{
    FnMated<-mapply(minCheck, round(colSums(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq))), mF)
    FnMated[is.na(FnMated)] = 0 
  }
  
  #B. Total number of births in each patch w density dependence, given avg litter size "k"
  sumLit<-(randpoisLitter(FnMated, k, FnMated))  #assumes min=2, max=9
  B<- sumLit*(1-((mXY + mXD + mXX)/CPop))
  B<- ifelse(B<0, 0, B)
  B <- (randpois(dsubpop = B, lambda = 1, subpop = B))
  
  #C. Mating fractions among all adult genotype pairs, assuming random mating within spatial patches
  allmate<- (mXX * (mXY + mXD))
  P1<- (mXX * mXY) /allmate
  P2<- (mXX * mXD) /allmate
  
  #D. calculate probability of offspring genotypes
  oXY<- P1*P1XY
  oXD<- P2*P2XD
  oXX<- P1*P1XX + P2*P2XX
  
  #E. Adjust offspring probabilities due to fitness costs. GD inds have reduced likelihood of birth
  if(length(PolyandrySeq)>1){
    oXD <- oXD-(oXD*(1-spD)*(sum(PolyandryProbs[-1])))
  }
  
  #F. Establish matrix of mating probabilities for randBirth2 function
  pBirth <-matrix(c(oXY,oXD,oXX),nrow=N*N,ncol=3,byrow=FALSE) #Male Bias. 
  
  #G. Establish final integer number of offpsring of each genotype
  Births <- randBirthXS(Birth = Birthnull, pBirth = pBirth, B = B)
  B1 <- Births[,1]
  B2 <- Births[,2]
  B3 <- Births[,3]
  
  
  ################################################################################
  #Step 3: Adult dispersal

  if((time.step-8)%%Dyr==0){
    
    M1 <- (rbinom(n = mXY, size=mXY, prob = dA))
    M2 <- (rbinom(n = mXD, size=mXD, prob = dA))
    M3 <- (rbinom(n = mXX, size=mXX, prob = dA))
    
    #Returns vector of new immigrants in each patch
    #Needs shape and rate parameters
    #Needs pSize parameter for informing range of potential destination patches
    tM1 <- (natMigration2(nMig=M1, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM2 <- (natMigration2(nMig=M2, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM3 <- (natMigration2(nMig=M3, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    
    #Update adult pops
    nXY <- mXY + tM1 - M1
    nXD <- mXD + tM2 - M2
    nXX <- mXX + tM3 - M3
    nXY <- ifelse(nXY<1, 0, nXY) #Ensure population does not go below 1. At least one needed to mate
    nXD <- ifelse(nXD<1, 0, nXD) 
    nXX <- ifelse(nXX<1, 0, nXX) 
  } else{
    nXY <- mXY
    nXD <- mXD
    nXX <- mXX
    nXY <- ifelse(nXY<1, 0, nXY) #Ensure population does not go below 1. At least one needed to mate
    nXD <- ifelse(nXD<1, 0, nXD) 
    nXX <- ifelse(nXX<1, 0, nXX) 
    
  }
  
  ##################################################################################
  #Step 4: Juvenile Dispersal
  
  if((time.step-8)%%Dyr==0){
    
    #Calc number of dispersers
    M4 <- (rbinom(n = B1, size=B1, prob = dJ))
    M5 <- (rbinom(n = B2, size=B2, prob = dJ))
    M6 <- (rbinom(n = B3, size=B3, prob = dJ))
    
    #Return number of new immigrants in each patch
    tM4 <- (natMigration2(nMig=M4, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM5 <- (natMigration2(nMig=M5, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    tM6 <- (natMigration2(nMig=M6, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
    
    #Update Juvenile populations from births and dispersal
    njXY <- jXY + B1 - M4 + tM4 
    njXD <- jXD + B2 - M5 + tM5
    njXX <- jXX + B3 - M6 + tM6
    
    njXY <- ifelse(njXY<0, 0,njXY)
    njXD <- ifelse(njXD<0, 0,njXD)
    njXX <- ifelse(njXX<0, 0,njXX)
    njXY[is.na(njXY)] = 0
    njXD[is.na(njXD)] = 0
    njXX[is.na(njXX)] = 0
  } else {
    njXY <- jXY + B1 
    njXD <- jXD + B2
    njXX <- jXX + B3
    
    njXY <- ifelse(njXY<0, 0,njXY)
    njXD <- ifelse(njXD<0, 0,njXD)
    njXX <- ifelse(njXX<0, 0,njXX)
    njXY[is.na(njXY)] = 0
    njXD[is.na(njXD)] = 0
    njXX[is.na(njXX)] = 0
  }
  ################################################################################
  #Step 5: Record important results
  
  #####
  #Sum births
  totBirth<-sum(B)
  
  #Step 5: Calculate end of time-step totals
  patchTotal <- nXY+nXD+nXX+njXY+njXD+njXX
  PopTot <- sum(patchTotal) #Global population
  
  GDMO <- nXD  #not including nRD because dead end individuals that wont drive themselves
  GDMOTot <- sum(GDMO)   #Global GDMO
  
  #Binary vector indicating if GDMOs are present (1) or not (0) in each patch
  Present <- ifelse(GDMO<1, 0,1)
  Present <- as.vector(Present)
  
  #Freq of gd allele in male pop
  GDMales<- (nXD)/(nXY+nXD)
  
  #genotype total for animation
  fqnXD<-sum(nXD+njXD)/PopTot
  fqnXY<-sum(nXY+njXY)/PopTot
  fqnXX<-sum(nXX+njXX)/PopTot
  
  #Sex Ratio in each patch (>0.5 is skewed male, <0.5 is skewed female)
  sexRatio<- (nXY+nXD+njXY+njXD) / patchTotal
  
  #Number of births in each patch
  #already contained in vector B
  
  #Out1: global frequency of gdmo individuals in context of total population
  #      vector of size NL*NL
  GDFrq <- GDMOTot/PopTot #Global Frequency
  GDFrq[is.na(GDFrq)] = 0

  #Out3: frequency of cells with at least 1 gdmo adult present (invaded)
  InvFrq <- sum(Present)/(N*N)
  
  #Out4: frequency of patches that have 0 individuals wt or gdmo (extinct)
  Extinct <- ifelse(patchTotal<1, 1, 0)
  ExtFrq <- sum(Extinct)/(N*N)
  
  
  #Tracking the mean,min,max distance from patches with wt present to nearest patch with gd
  wtTotal<-nXY+nXX+njXY+njXX
  wtPresent <- ifelse(wtTotal<1, 0,1)
  wt2gdDists<-dist2GdPatch(wtPresent, Present, N)
  
  #Tracking the time step...maybe not necessary
  ntime.step <-(time.step+1)
  
  dx <- c(nXY,nXD,nXX,njXY,njXD,njXX,
          patchTotal, sexRatio, B, GDMales,
          totBirth, GDMOTot, GDFrq, PopTot,
          InvFrq, ExtFrq, fqnXD, fqnXY, fqnXX, wt2gdDists, ntime.step) #vector of populations.
  
  #######################################
  ##combine results into a single vector dx
  ######################################
  list(dx)  
}

#Y-shredder gene drive dynamics. Juvenile dispersal frequency mirrors adult dispersal frequency parameter. Tracks chasing dynamics
Yshred.v2.chase <- function(t, y, parms, signal){
  
  ################################################################################  
  #Setup: Load Parameters and initital abundances
  
  N<-NL
  g=parms[[1]]            #Development rate
  u=parms[[2]]            #Mortality rate
  k=parms[[3]]            #Litter size
  v=parms[[4]]            #Female mating rate
  CPop=parms[[5]]         #Carrying capacity of each patch
  Andry=parms[[6]]        #Max male mates per female
  Gyny=parms[[7]]         #Max female mates per male
  kerShape=parms[[8]]     #Dispersal kernel shape
  kerRate=parms[[9]]      #dispersal kernel rate
  dJ=parms[[10]]          #Proportion of juveniles dispersing
  dA=parms[[10]]          #Proportion of adults dispersing
  gdInt=parms[[12]]       #number of gene drive individuals introduced
  gdSupp=parms[[13]]      #release schedule (every 1-10yrs), or once
  intLoc=parms[[14]]      #release all in middle (1) or each ind randomly (2)
  Pc=parms[[15]]          #prob of HEG1 cutting at homologous locus
  Pn=parms[[16]]          #prob of NHEJ after HEG1 cut (i.e. successful conversion)
  Py=parms[[17]]          #prob of HEG2 cutting at Y chromosome
  Xpass=parms[[18]]       #prob of XO female passing X chrom
  s=parms[[19]]           #Fitness cost for homozygous GDMO
  sh=parms[[20]]          #Fitness cost for heterozygous GDMO
  
  #Setup: Create vector of N from vector loaded into model
  XXw <- y[(0*(N^2) + 1): (1*(N^2))] 
  XYw <- y[(1*(N^2) + 1): (2*(N^2))]
  XOw <- y[(2*(N^2) + 1): (3*(N^2))]
  XXH <- y[(3*(N^2) + 1): (4*(N^2))]
  XXh <- y[(4*(N^2) + 1): (5*(N^2))]
  XYH <- y[(5*(N^2) + 1): (6*(N^2))]
  XYh <- y[(6*(N^2) + 1): (7*(N^2))]
  XOH <- y[(7*(N^2) + 1): (8*(N^2))]
  XOh <- y[(8*(N^2) + 1): (9*(N^2))]
  jXXw <- y[(9*(N^2) + 1): (10*(N^2))]
  jXYw <- y[(10*(N^2) + 1): (11*(N^2))]
  jXOw <- y[(11*(N^2) + 1): (12*(N^2))]
  jXXH <- y[(12*(N^2) + 1): (13*(N^2))]
  jXXh <- y[(13*(N^2) + 1): (14*(N^2))]
  jXYH <- y[(14*(N^2) + 1): (15*(N^2))]
  jXYh <- y[(15*(N^2) + 1): (16*(N^2))]
  jXOH <- y[(16*(N^2) + 1): (17*(N^2))]
  jXOh <- y[(17*(N^2) + 1): (18*(N^2))]
  time.step<-y[((22*(N^2))+23)]
  
  ################################################################################  
  #GD Release at specific time steps
  
  #Scheduled release of GD individuals (heterozygous females), once here, successively below
  if(gdSupp==11){
    if(time.step==7){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- gdInt #Add specific n of GD ind to central patch 
      }
      if(intLoc==2){              #Random placement
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
  }
  #First release when supplemental releases scheduled
  if(gdSupp<11) {
    if(time.step==7){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- gdInt #Always add at the beginning
      }
      if(intLoc==2){              #Random placement
        for(i in 1:gdInt){
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
    #Supplemental release every X years
    if(time.step>7 & (((time.step-7)/6)%%gdSupp)==0){
      if(intLoc==1){
        XXH[(N^2/2)+0.5] <- (XXH[(N^2/2)+0.5]+gdInt) #Add new GD individuals to central patch
      }
      if(intLoc==2){
        for(i in 1:gdInt){              #Random placement
          loc<-sample(N^2,1,)
          XXH[loc]<- (XXH[loc]+1)
        }
      }
    }
  }
  
  ################################################################################
  #Step 1: Mortality and Maturation
  
  #Determine number of adult mortalities
  D1 <- (rbinom(n = XXw, prob = u, size = XXw))
  D2 <- (rbinom(n = XYw, prob = u, size = XYw))
  D3 <- (rbinom(n = XXH, prob = u, size = XXH))
  D4 <- (rbinom(n = XXh, prob = u, size = XXh))
  D5 <- (rbinom(n = XYH, prob = u, size = XYH))
  D6 <- (rbinom(n = XYh, prob = u, size = XYh))
  D7 <- (rbinom(n = XOH, prob = u, size = XOH))
  D8 <- (rbinom(n = XOh, prob = u, size = XOh))
  D9 <- (rbinom(n = XOw, prob = u, size = XOw))
  
  #Determine number of juvenile mortalities
  D10 <- (rbinom(n = jXXw, prob = u, size = jXXw))
  D11 <- (rbinom(n = jXYw, prob = u, size = jXYw))
  D12 <- (rbinom(n = jXXH, prob = u, size = jXXH))
  D13 <- (rbinom(n = jXXh, prob = u, size = jXXh))
  D14 <- (rbinom(n = jXYH, prob = u, size = jXYH))
  D15 <- (rbinom(n = jXYh, prob = u, size = jXYh))
  D16 <- (rbinom(n = jXOH, prob = u, size = jXOH))
  D17 <- (rbinom(n = jXOh, prob = u, size = jXOh))
  D18 <- (rbinom(n = jXOw, prob = u, size = jXOw))
  
  #Update juvenile population after mortality
  jXXw <- jXXw - D10
  jXYw <- jXYw - D11
  jXXH <- jXXH - D12
  jXXh <- jXXh - D13
  jXYH <- jXYH - D14
  jXYh <- jXYh - D15
  jXOH <- jXOH - D16
  jXOh <- jXOh - D17
  jXOw <- jXOw - D18
  jXXw <- ifelse(jXXw<0,0,jXXw) 
  jXYw <- ifelse(jXYw<0,0,jXYw) 
  jXXH <- ifelse(jXXH<0,0,jXXH) 
  jXXh <- ifelse(jXXh<0,0,jXXh) 
  jXYH <- ifelse(jXYH<0,0,jXYH) 
  jXYh <- ifelse(jXYh<0,0,jXYh) 
  jXOH <- ifelse(jXOH<0,0,jXOH) 
  jXOh <- ifelse(jXOh<0,0,jXOh) 
  jXOw <- ifelse(jXOw<0,0,jXOw) 
  jXXw[is.na(jXXw)] = 0
  jXYw[is.na(jXYw)] = 0
  jXXH[is.na(jXXH)] = 0
  jXXh[is.na(jXXh)] = 0
  jXYH[is.na(jXYH)] = 0
  jXYh[is.na(jXYh)] = 0
  jXOH[is.na(jXOH)] = 0
  jXOh[is.na(jXOh)] = 0
  jXOw[is.na(jXOw)] = 0
  
  #Poisson of maturation
  J1 <- (rbinom(n = jXXw, prob = g, size = jXXw))
  J2 <- (rbinom(n = jXYw, prob = g, size = jXYw))
  J3 <- (rbinom(n = jXXH, prob = g, size = jXXH))
  J4 <- (rbinom(n = jXXh, prob = g, size = jXXh))
  J5 <- (rbinom(n = jXYH, prob = g, size = jXYH))
  J6 <- (rbinom(n = jXYh, prob = g, size = jXYh))
  J7 <- (rbinom(n = jXOH, prob = g, size = jXOH))
  J8 <- (rbinom(n = jXOh, prob = g, size = jXOh))
  J9 <- (rbinom(n = jXOw, prob = g, size = jXOw))
  
  #Update adult subpopulations
  mXXw <- XXw + J1 - D1
  mXYw <- XYw + J2 - D2
  mXXH <- XXH + J3 - D3
  mXXh <- XXh + J4 - D4
  mXYH <- XYH + J5 - D5
  mXYh <- XYh + J6 - D6
  mXOH <- XOH + J7 - D7
  mXOh <- XOh + J8 - D8
  mXOw <- XOw + J9 - D9
  mXXw <- ifelse(mXXw<0,0,mXXw) 
  mXYw <- ifelse(mXYw<0,0,mXYw) 
  mXXH <- ifelse(mXXH<0,0,mXXH) 
  mXXh <- ifelse(mXXh<0,0,mXXh) 
  mXYH <- ifelse(mXYH<0,0,mXYH) 
  mXYh <- ifelse(mXYh<0,0,mXYh)
  mXOH <- ifelse(mXOH<0,0,mXOH)
  mXOh <- ifelse(mXOh<0,0,mXOh)
  mXOw <- ifelse(mXOw<0,0,mXOw)
  mXXw[is.na(mXXw)] = 0
  mXYw[is.na(mXYw)] = 0
  mXXH[is.na(mXXH)] = 0
  mXXh[is.na(mXXh)] = 0
  mXYH[is.na(mXYH)] = 0
  mXYh[is.na(mXYh)] = 0  
  mXOH[is.na(mXOH)] = 0 
  mXOh[is.na(mXOh)] = 0 
  mXOw[is.na(mXOw)] = 0 
  
  #Update juvenile population after maturation
  jXXw <- jXXw - J1
  jXYw <- jXYw - J2
  jXXH <- jXXH - J3
  jXXh <- jXXh - J4
  jXYH <- jXYH - J5
  jXYh <- jXYh - J6
  jXOH <- jXOH - J7
  jXOh <- jXOh - J8
  jXOw <- jXOw - J9
  jXXw <- ifelse(jXXw<0,0,jXXw) 
  jXYw <- ifelse(jXYw<0,0,jXYw) 
  jXXH <- ifelse(jXXH<0,0,jXXH) 
  jXXh <- ifelse(jXXh<0,0,jXXh) 
  jXYH <- ifelse(jXYH<0,0,jXYH) 
  jXYh <- ifelse(jXYh<0,0,jXYh) 
  jXOH <- ifelse(jXOH<0,0,jXOH) 
  jXOh <- ifelse(jXOh<0,0,jXOh) 
  jXOw <- ifelse(jXOw<0,0,jXOw) 
  jXXw[is.na(jXXw)] = 0
  jXYw[is.na(jXYw)] = 0
  jXXH[is.na(jXXH)] = 0
  jXXh[is.na(jXXh)] = 0
  jXYH[is.na(jXYH)] = 0
  jXYh[is.na(jXYh)] = 0
  jXOH[is.na(jXOH)] = 0
  jXOh[is.na(jXOh)] = 0
  jXOw[is.na(jXOw)] = 0
  
  ################################################################################
  #Step 2: Adult mating
  
  #A. Mating population established (mating males "mM" and mating females "mF")
  mM <- mXYw + mXYH + mXYh                    #Male mating population
  mF <- v*(mXXw + mXXH + mXXh + mXOH + mXOh + mXOw)  #Female Mating Population 
  
  # Determine number mated females given max number of mates per males and females
  #Requires polyandry probabilities (likelihood of n mates for females for 1:max)
  if(Andry==1){
    FnMated<-mapply(minCheck, round(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq)), mF)
    FnMated[is.na(FnMated)] = 0 
  }else{
    FnMated<-mapply(minCheck, round(colSums(sapply(mM, MatedFemByPoly, MaxMatesPerMale=Gyny, PolyProbs=PolyandryProbs, PolySeq=PolyandrySeq))), mF)
    FnMated[is.na(FnMated)] = 0 
  }
  
  #B. Total number of births in each patch w density dependence, given avg litter size "k"
  sumLit<-(randpoisLitter(FnMated, k, FnMated))  #assumes min=2, max=9
  B<- sumLit*(1-((mXXw + mXYw + mXXH + mXXh + mXYH + mXYh + mXOH + mXOh + mXOw)/CPop))
  B<- ifelse(B<0, 0,B)
  B <- (randpois(dsubpop = B, lambda = 1, subpop = B))

  
  #C. Mating fractions among all adult genotype pairs, assuming random mating within spatial patches
  allmate<- (mXXw * (mXYw + mXYH +mXYh)) + (mXXH * (mXYw + mXYH +mXYh)) + (mXXh * (mXYw + mXYH +mXYh)) + (mXOH * (mXYw + mXYH +mXYh)) + (mXOh * (mXYw + mXYH +mXYh)) + (mXOw * (mXYw + mXYH +mXYh))
  P1<-(mXXw * mXYw) / allmate
  P2<-(mXXw * mXYh) / allmate
  P3<-(mXXw * mXYH) / allmate
  P4<-(mXXh * mXYw) / allmate
  P5<-(mXXh * mXYh) / allmate
  P6<-(mXXh * mXYH) / allmate
  P7<-(mXXH * mXYw) / allmate
  P8<-(mXXH * mXYh) / allmate
  P9<-(mXXH * mXYH) / allmate
  P10<-(mXOh * mXYw) / allmate
  P11<-(mXOh * mXYh) / allmate
  P12<-(mXOh * mXYH) / allmate
  P13<-(mXOH * mXYw) / allmate
  P14<-(mXOH * mXYh) / allmate
  P15<-(mXOH * mXYH) / allmate
  P16<-(mXOw * mXYw) / allmate
  P17<-(mXOw * mXYH) / allmate
  P18<-(mXOw * mXYh) / allmate
  
  #D. Offspring probs sum of mating fraction and inheritance calculation for all possible parent combos that could produce offspring geno
  oXXw<- P1*P1XXw + P2*P2XXw + P4*P4XXw + P5*P5XXw + P10*P10XXw + P11*P11XXw + P15*P15XXw + P18*P18XXw
  oXYw<- P1*P1XYw + P2*P2XYw + P4*P4XYw + P5*P5XYw + P10*P10XYw + P11*P11XYw + P15*P15XYw + P18*P18XYw
  oXXH<- P2*P2XXH + P3*P3XXH + P4*P4XXH + P5*P5XXH + P6*P6XXH + P7*P7XXH + P8*P8XXH + P9*P9XXH + P10*P10XXH + P11*P11XXH + P12*P12XXH + P13*P13XXH + P14*P14XXH + P15*P15XXH + P17*P17XXH + P18*P18XXH
  oXXh<- P2*P2XXh + P3*P3XXh + P4*P4XXh + P5*P5XXh + P6*P6XXh + P7*P7XXh + P8*P8XXh + P10*P10XXh + P11*P11XXh + P12*P12XXh + P13*P13XXh + P14*P14XXh + P17*P17XXh + P18*P18XXh
  oXYH<- P2*P2XYH + P3*P3XYH + P4*P4XYH + P5*P5XYH + P6*P6XYH + P7*P7XYH + P8*P8XYH + P9*P9XYH + P10*P10XYH + P11*P11XYH + P12*P12XYH + P13*P13XYH + P14*P14XYH + P15*P15XYH + P17*P17XYH + P18*P18XYH
  oXYh<- P2*P2XYh + P3*P3XYh + P4*P4XYh + P5*P5XYh + P6*P6XYh + P7*P7XYh + P8*P8XYh + P10*P10XYh + P11*P11XYh + P12*P12XYh + P13*P13XYh + P14*P14XYh + P17*P17XYh + P18*P18XYh
  oXOH<- P2*P2XOH + P3*P3XOH + P4*P4XOH + P5*P5XOH + P6*P6XOH + P7*P7XOH + P8*P8XOH + P9*P9XOH + P10*P10XOH + P11*P11XOH + P12*P12XOH + P13*P13XOH + P14*P14XOH + P15*P15XOH + P17*P17XOH + P18*P18XOH
  oXOh<- P2*P2XOh + P3*P3XOh + P4*P4XOh + P5*P5XOh + P6*P6XOh + P7*P7XOh + P8*P8XOh + P10*P10XOh + P11*P11XOh + P12*P12XOh + P13*P13XOh + P14*P14XOh + P17*P17XOh + P18*P18XOh
  oXOw<- P10*P10XOw + P11*P11XOw + P16*P16XOw
  #all OY and failed conversions are nonviable
  onv <- P2*P2nv + P3*P3nv + P4*P4nv + P5*P5nv + P6*P6nv + P7*P7nv + P8*P8nv + P10*P10nv + P11*P11nv + P12*P12nv + P13*P13nv + P14*P14nv + P17*P17nv + P18*P18nv
  
  #E. Establish matrix of birth probabilities
  pBirth<-matrix(c(oXXw,oXYw,oXXH,oXXh,oXYH,oXYh,oXOH,oXOh,oXOw,onv), nrow=N*N, ncol=9, byrow=F)
  pBirth[is.na(pBirth)] = 0
  
  #F. Adjust offspring probabilities due to fitness costs. GD inds have reduced likelihood of birth
  oXXH<- oXXH*(1-s)
  oXYH<- oXYH*(1-s)
  oXOH<- oXOH*(1-s)
  oXXh<- oXXh*(1-sh)
  oXYh<- oXYh*(1-sh)
  oXOh<- oXOh*(1-sh)
  
  #G. Establish final integer number of offspring of each genotype
  Births <- randBirthYS(Birth = Birthnull, pBirth = pBirth, B = B)
  B1 <- Births[,1]
  B2 <- Births[,2]
  B3 <- Births[,3]
  B4 <- Births[,4]
  B5 <- Births[,5]
  B6 <- Births[,6]
  B7 <- Births[,7]
  B8 <- Births[,8]
  B9 <- Births[,9]
  
  ################################################################################
  #Step 3: Adult dispersal
  
  #Determine number of dispersers
  M1 <- (rbinom(n = mXXw, prob = dA, size = mXXw))
  M2 <- (rbinom(n = mXYw, prob = dA, size = mXYw))
  M3 <- (rbinom(n = mXXH, prob = dA, size = mXXH))
  M4 <- (rbinom(n = mXXh, prob = dA, size = mXXh))
  M5 <- (rbinom(n = mXYH, prob = dA, size = mXYH))
  M6 <- (rbinom(n = mXYh, prob = dA, size = mXYh))
  M7 <- (rbinom(n = mXOH, prob = dA, size = mXOH))
  M8 <- (rbinom(n = mXOh, prob = dA, size = mXOh))
  M9 <- (rbinom(n = mXOw, prob = dA, size = mXOw))
  
  #Returns vector of new immigrants in each patch
  tM1 <- (natMigration2(nMig=M1, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM2 <- (natMigration2(nMig=M2, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM3 <- (natMigration2(nMig=M3, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM4 <- (natMigration2(nMig=M4, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM5<-  (natMigration2(nMig=M5, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM6 <- (natMigration2(nMig=M6, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))   
  tM7 <- (natMigration2(nMig=M7, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM8 <- (natMigration2(nMig=M8, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM9 <- (natMigration2(nMig=M9, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  
  #Update adult pops
  nXXw <- mXXw + tM1 - M1
  nXYw <- mXYw + tM2 - M2
  nXXH <- mXXH + tM3 - M3
  nXXh <- mXXh + tM4 - M4
  nXYH <- mXYH + tM5 - M5
  nXYh <- mXYh + tM6 - M6
  nXOH <- mXOH + tM7 - M7
  nXOh <- mXOh + tM8 - M8
  nXOw <- mXOw + tM9 - M9
  #Ensure population does not go below 1. At least one needed to mate
  nXXw <- ifelse(nXXw<1, 0, nXXw) 
  nXYw <- ifelse(nXYw<1, 0, nXYw) 
  nXXH <- ifelse(nXXH<1, 0, nXXH) 
  nXXh <- ifelse(nXXh<1, 0, nXXh) 
  nXYH <- ifelse(nXYH<1, 0, nXYH) 
  nXYh <- ifelse(nXYh<1, 0, nXYh) 
  nXOH <- ifelse(nXOH<1, 0, nXOH) 
  nXOh <- ifelse(nXOh<1, 0, nXOh)
  nXOw <- ifelse(nXOw<1, 0, nXOw)
  
  ################################################################################
  #Step 4: Juvenile dispersal
  
  #Calc number of dispersers
  M10 <- (rbinom(n = B1, prob = dJ, size = B1))
  M11 <- (rbinom(n = B2, prob = dJ, size = B2))
  M12 <- (rbinom(n = B3, prob = dJ, size = B3))
  M13 <- (rbinom(n = B4, prob = dJ, size = B4))
  M14 <- (rbinom(n = B5, prob = dJ, size = B5))
  M15 <- (rbinom(n = B6, prob = dJ, size = B6))
  M16 <- (rbinom(n = B7, prob = dJ, size = B7))
  M17 <- (rbinom(n = B8, prob = dJ, size = B8))
  M18 <- (rbinom(n = B9, prob = dJ, size = B9))
  
  #Return number of new immigrants in each patch
  tM10 <- (natMigration2(nMig=M10, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM11 <- (natMigration2(nMig=M11, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM12 <- (natMigration2(nMig=M12, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM13 <- (natMigration2(nMig=M13, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate)) 
  tM14 <- (natMigration2(nMig=M14, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM15 <- (natMigration2(nMig=M15, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM16 <- (natMigration2(nMig=M16, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM17 <- (natMigration2(nMig=M17, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  tM18 <- (natMigration2(nMig=M18, N=NL, pSize=0.5, kernShape=kerShape, kernRate=kerRate))
  
  #Update Juvenile populations prior to mortality/maturation
  njXXw <- jXXw + B1 - M10 + tM10
  njXYw <- jXYw + B2 - M11 + tM11
  njXXH <- jXXH + B3 - M12 + tM12
  njXXh <- jXXh + B4 - M13 + tM13
  njXYH <- jXYH + B5 - M14 + tM14
  njXYh <- jXYh + B6 - M15 + tM15
  njXOH <- jXOH + B7 - M16 + tM16
  njXOh <- jXOh + B8 - M17 + tM17
  njXOw <- jXOw + B9 - M18 + tM18
  
  #Setting negative values and NAs to zero
  njXXw <- ifelse(njXXw<0, 0,njXXw)
  njXYw <- ifelse(njXYw<0, 0,njXYw)
  njXXH <- ifelse(njXXH<0, 0,njXXH)
  njXXh <- ifelse(njXXh<0, 0,njXXh)
  njXYH <- ifelse(njXYH<0, 0,njXYH)
  njXYh <- ifelse(njXYh<0, 0,njXYh)
  njXOH <- ifelse(njXOH<0, 0,njXOH)
  njXOh <- ifelse(njXOh<0, 0,njXOh)
  njXOw <- ifelse(njXOw<0, 0,njXOw)
  njXXw[is.na(njXXw)] = 0
  njXYw[is.na(njXYw)] = 0
  njXXH[is.na(njXXH)] = 0
  njXXh[is.na(njXXh)] = 0
  njXYH[is.na(njXYH)] = 0
  njXYh[is.na(njXYh)] = 0
  njXOH[is.na(njXOH)] = 0
  njXOh[is.na(njXOh)] = 0
  njXOw[is.na(njXOw)] = 0
  
  ################################################################################
  #Step 5: Record important results
  
  #Sum births
  totBirth<-sum(B)
  
  #Step 5: Calculate end of time-step totals
  subTotal <- nXXw+nXYw+nXOw+nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXw+njXYw+njXOw+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh
  PopTot <- sum(subTotal) #Global population
  
  GDMO <- nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh
  GDMOTot <- sum(GDMO)   #Global GDMO
  
  #Binary vector indicating if GDMOs are present (1) or not (0) in each patch
  Present <- ifelse(GDMO<1,0,1)
  Present <- as.vector(Present)
  
  #Freq of gd allele in male pop
  GDPatch<- (nXXH+nXXh+nXYH+nXYh+nXOH+nXOh+njXXH+njXXh+njXYH+njXYh+njXOH+njXOh)/subTotal
  
  #Sex Ratio
  sexRatio<- (nXYw+nXYH+nXYh+njXYw+njXYH+njXYh)/subTotal
  
  #Genotype frequencies
  wtFq<-sum(nXXw+nXYw+nXOw+njXXw+njXYw+njXOw)/PopTot
  hetFq<-sum(nXXh+nXYh+nXOh+njXXh+njXYh+njXOh)/PopTot
  homFq<-sum(nXXH+nXYH+nXOH+njXXH+njXYH+njXOH)/PopTot
  
  #Specific genotype frequencies
  XXwFq<-sum(nXXw + njXXw)/PopTot
  XYwFq<-sum(nXYw + njXYw)/PopTot
  XOwFq<-sum(nXOw + njXOw)/PopTot
  XXHFq<-sum(nXXH + njXXH)/PopTot
  XYHFq<-sum(nXYH + njXYH)/PopTot
  XOHFq<-sum(nXOH + njXOH)/PopTot
  XXhFq<-sum(nXXh + njXXh)/PopTot
  XYhFq<-sum(nXYh + njXYh)/PopTot
  XOhFq<-sum(nXOh + njXOh)/PopTot
  
  
  #Out1: global frequency of gdmo individuals in context of total population
  #      vector of size NL*NL
  GDFrq <- GDMOTot/PopTot #Global Frequency
  GDFrq[is.na(GDFrq)] = 0
  
  #Out2: resistant individual freq
  ResFrq<-0
  
  #Out2: Global abundance over time, single value per time step
  # just use "PopTot"...don't need to make new object
  
  #Out3: frequency of cells with at least 1 gdmo adult present (invaded)
  InvFrq <- sum(Present)/(N*N)
  
  #Out4: frequency of patches that have 0 individuals wt or gdmo (extinct)
  Extinct <- ifelse(subTotal<1, 1, 0)
  ExtFrq <- sum(Extinct)/(N*N)
  
  #True Female Mating Rate
  TrueFMR<-FnMated/(mXXw + mXXH + mXXh + mXOH + mXOh + mXOw)
  
  #Tracking the mean,min,max distance from patches with wt present to nearest patch with gd
  wtTotal<-nXXw+nXYw+nXOw+njXXw+njXYw+njXOw
  wtPresent <- ifelse(wtTotal<1, 0,1)
  wt2gdDists<-dist2GdPatch(wtPresent, Present, N)
  
  
  #Tracking the time step...maybe not necessary
  ntime.step <-(time.step+1)
  
  dx <- c(nXXw,nXYw,nXOw,nXXH,nXXh,nXYH,nXYh,nXOH,nXOh,njXXw,njXYw,njXOw,njXXH,njXXh,njXYH,njXYh,njXOH,njXOh,
          subTotal,sexRatio, B, GDPatch,
          totBirth,GDMOTot,GDFrq,ResFrq,PopTot,
          InvFrq,ExtFrq,wt2gdDists,wtFq,hetFq,homFq,XXwFq,XYwFq,XOwFq,XXHFq,XYHFq,XOHFq,XXhFq,XYhFq,XOhFq,ntime.step) #vector of populations. 
  
  ##combine results into a single vector dx
  list(dx)    
}

#Determines XY position of patch centroids within square patch lattice.
calculateCentroids<-function(n)
{
  xRow <- seq(from = 0.5, to = (n - 0.5), by = 1)
  X <- matrix(rep(xRow, n), nrow = n, byrow = T)
  X <- as.vector(t(X)) # flatten out matrix, row by row
  
  yCol <- seq(from = (n - 0.5), to = 0.5, by = -1)
  Y <- matrix(rep(yCol, n), nrow = n, byrow = F)
  Y <- as.vector(t(Y)) # flatten out matrix, row by row
  
  ctr<-matrix(c(X,Y), nrow=(n*n), byrow=F)
  
  return(ctr)
}

#Determines final abundance in each patch following migration
natMigration2<-function(nMig, N , pSize, kernShape, kernRate) 
{
  cent<-calculateCentroids(N) #calculate the centroid locations based on lattice side length
  pmove<-c()
  tmove<-c()
  for (j in 1:length(nMig)) 
  {
    wMig<-nMig[j]
    if(wMig>0){
      patDists<- sqrt( (cent[,1] - (cent[j,1]))^2 + (cent[,2]- (cent[j,2]))^2) #euclidean distance between patch j and all other centroids
      migDists<-rgamma(wMig, shape=kernShape, rate=kernRate) #sample dispersal kernel for vector of migration distances
      migDists[migDists<1]<-1
      migDists[migDists>(max(patDists))]<-max(patDists)
      for(i in 1:length(migDists))
      {
        diff<- abs(patDists-migDists[i]) #diff between distance to patches and migration distance of each disperser
        opts<-which(diff < pSize & diff > 0) #list of possible patches, those within specified distance threshold
        if(length(opts)>1){ #when there are more than 1 patch options
          choice<-sample(opts, 1) #sample the options and choose one at random
          pmove<-c(pmove, choice) #add destination patch to the vector of migrant destinations
        }
        else if(length(opts)==1){  #if there is only one option
          pmove<-c(pmove,opts)       #use that single option and add as destination patch
        }
        else if(length(opts)==0){ #
          pmove<-c(pmove,j)
        }
      }
    }
  }
  
  for (k in 1:length(nMig))
  {
    num<-length(which(pmove == k))
    tmove<-c(tmove, num)
  }
  return(tmove) 
}

#Determines discrete number of births of each potential genotype in the Y-Shredder system using weighted sampling.
randBirthYS <- function(Birth, pBirth, B)
{
  gens<-length(pBirth[1,])
  for (i in 1:nrow(Birth))
  {
    if(B[i] > 0){
      prob <- pBirth[i,]
      nsample <- B[i]
      
      #random sample
      sampidx <- sample(gens,
                        nsample,
                        prob=prob,
                        replace=TRUE)
      
      freqTable <- tabulate(sampidx, nbins=10)
      Birth[i,] <- freqTable
    }
  }
  return (Birth)
}

##Determines discrete number of births of each potential genotype in the X-Shredder system using weighted sampling.
randBirthXS <- function(Birth, pBirth, B)
{
  gens<-length(pBirth[1,])
  for (i in 1:nrow(Birth))
  {
    if(B[i] > 0){
      prob <- pBirth[i,]
      nsample <- B[i]
      
      #random sample
      sampidx <- sample(gens,
                        nsample,
                        prob=prob,
                        replace=TRUE)
      
      freqTable <- tabulate(sampidx, nbins=3)
      Birth[i,] <- freqTable
    }
  }
  return (Birth)
}

#Determines discrete numbers of individuals experiencing dispersal or demographic processes for each patch, given relevant parameter and patch abundances.
randpois <- function(dsubpop,lambda,subpop)
{
  for (i in 1:length(dsubpop))
  {
    lambda <- lambda #
    #If statement in case sample vector is 0, this way data.table will work. Below added values are replaced with 0
    tsample <- subpop[i]
    nsample <- ifelse(tsample<1, 1,tsample)
    #random sample rpois(trials,lambda)
    sampidx <- rpois(nsample,lambda)
    count <- sum(sampidx)
    
    #If statement in case subpopulation was 0
    Refadj <- ifelse(tsample==0, 0,1)
    count <- count*Refadj
    count <- as.numeric(count)
    
    #Update by referencing cells in vector (Move) and adding Ref value
    dsubpop[i] <- count
  }
  return (dsubpop)
}

#Determines discrete number of individuals in each litter, capped at minimum and maximum values for house mice
randpoisLitter <- function(dsubpop,lambda,subpop)
{
  for (i in 1:length(dsubpop))
  {
    lambda <- lambda #
    #If statement in case sample vector is 0, this way data.table will work. Below added values are replaced with 0
    tsample <- subpop[i]
    nsample <- ifelse(tsample<1, 1,tsample)
    #random sample rpois(trials,lambda)
    sampidx <- extraDistr::rtpois(nsample,lambda, 1, 9)
    count <- sum(sampidx)
    
    #If statement in case subpopulation was 0
    Refadj <- ifelse(tsample==0, 0,1)
    count <- count*Refadj
    count <- as.numeric(count)
    
    #Update by referencing cells in vector (Move) and adding Ref value
    dsubpop[i] <- count
  }
  return (dsubpop)
}

#Determines the number of mated females in each patch given male availability and mating system characteristics.
MatedFemByPoly<-function(x,MaxMatesPerMale,PolyProbs,PolySeq){
  ((x*MaxMatesPerMale)*PolyandryProbs)/PolySeq
}

#Determines minimum value
minCheck<-function(x,y){
  min(x,y)
}

#Calculates the distance between wt individuals and nearest patch containing gene drive
dist2GdPatch<-function(wtPres, gdPres, N){
  if(sum(wtPres)==0){
    return(c(NA,NA,NA))
  } else if(sum(gdPres)==0) {
    return(c(NA,NA,NA))
  } else
    minDists<-c()
  cent<-calculateCentroids(N)
  gdPresLogic<-as.logical(gdPres)
  for(i in 1:length(wtPres)){
    if(wtPres[i] == 1){
      patDists<- sqrt( (cent[,1] - (cent[i,1]))^2 + (cent[,2]- (cent[i,2]))^2) #distance to all patches from wt-pres patch
      gdDists<-patDists[gdPresLogic] #distance to any patch with gd-pres
      minDists<-c(minDists,min(gdDists)) #collect the min distance to a gd-pres patch
    }
  }
  meanGdDist<-mean(minDists) # record average dist b/w wt-pres and nearest gd-pres
  minGdDist<-min(minDists) #record smallest dist b/w wt-pres and nearest gd-pres
  maxGdDist<-max(minDists) #record largest distance b/w wt-pres and nearest gd-pres
  distRes<-c(meanGdDist, minGdDist, maxGdDist)
  return(distRes)
}