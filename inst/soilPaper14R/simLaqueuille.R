isPaperBGC <- ("paperBGC" %in% commandArgs(trailingOnly = TRUE))
#isPaperBGC <- TRUE

library(ggplot2)
library(grid)   #unit
library(reshape)  # melt
library(RColorBrewer) # brewer.pal
library(twDEMC)
library(twDEMCPlot)
library(FME)

# simulating CO2 increase and bare soil, based on modEezy5
# based on 
if( isPaperBGC){
    library(twDev)
    loadPkg()
    baseFontSize <- 9  # pubs
    baseLineSize <- 0.6
	themeDefault <- theme_classic(base_size=baseFontSize)
} else {
    baseFontSize <- 16  # presentations
    baseLineSize <- 1.2
	themeDefault <- theme_classic(base_size=baseFontSize)
}

# gC/m2 and gN/m2, /yr
# From Perveen 2014 Symphony
B00 <- 16
A_Perveen14 <- 0.0317917*365 # decomposer cosumption rate of SOM (table 1, eq. 
L0 <- L0_Perveen14 <- 635
SOM0_Allard07 <- 12000
R0 <-  SOM0_Allard07 - L0_Perveen14   #640*0.1*0.18*1000 rho*10cm*Corg

NLoss_A_Perveen14 <- 0.0055 # phi_i - l*N - phi_up


parms0 <- list(
        cnB = 11         ##<< Perveen14, slightly higher than cnR
        ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        ,cnIR = 11-0.1    ##<< between micr and enzyme signal 
        #,cnIR = 4.5    ##<< between micr and enzyme signal 
        ,cnIL = 70      ##<< N poor substrate, Perveen14
        #,kN= 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
        #,kN= 1/(1/12)   ##<< enzyme tvr 1 month (Blagodatskaya 60 days priming)
        ,kN = 60       ##<< /yr enzyme turnover 60 times a year, each 6 days -> fast priming pulses
        ,km = 0.05       ##<< enzyme half-saturation constant, in magnitude of enzymes, determined by kN
        ,kNB = 0.8      ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R)
        ,kR = 1/(10)    ##<< 1/(x years) 
        #,kR = A_Perveen14 * B00 / (R0 *0.5)  #1/(10)     ##<< 1/(x years)     estimated from A in Perveen14
        #,kL = 2.4 #1/(0.33)   ##<< 1/(x years)     # formerly 1 year  # 4.23e-4*365
        ,kL = 2 # 1/yr, turns over X per year 
        #,aE = 1e-4*365 #0.001*365  ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,aE = 1e-3*365   ##<< C biomass allocated to enzymes gC/yr/microbial biomass, corresponds about 2% of uptake, toCheck
        #,km = 0.08       ##<< enzyme half-saturation constant
        #,km = 1           ##<< enzyme half-saturation constant, large value to make sure in linear part
        ,m = 0 # 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = (0.016906/0.8)*365 #1/60*365  ##<< biomass turnover rate from Perveen (s), but accounting for epsTvr tauPerv = epsTvr*tau
        #,tau = 1/10*365 #      ##<< biomass turnover rate of x days, with too low turnover time, too low microbial biomass
        #,tau = (0.016906+0.0368857)*365 #1/60*365  ##<< biomass turnover rate from Perveen (s+r) to match similar microbes pool
        #,tau = 1/20*365 #      ##<< biomass turnover rate of x days 
        ,eps = 0.33        ##<< carbon use efficiency  Perveen14: r/(s+r) 31% (here need to addionally account for 2% for enzymes)
        #,epsTvr = 1 #0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
        ,epsTvr =  0.8   ##<< carbon use efficiency of microbial tvr (predators respire)
        ,iR = 0         ##<< input by feedback modelled explicitely
        ,iL = 0.00505757*365*525         # g/m2 input per year mp*Cp
        ,plantNUp = 0                    # rate of organic N uptake
        ,useFixedAlloc=FALSE    ##<< set to true to use Fixed enzyme allocation (alpha = 0.5)
        ,nu=0.7                 ##<< aggregated microbial organic nitrogen assimilation efficiency (Manzoni 2008)
        #,kIP = 10.57 #0.0289652*365          ##<< plant uptake kIP I [1/yr]        
        #,iB = 0.0110068*365   ##<< immobilization flux [1/yr], from Perveen14 Table 1(i), is larger than required immobilization flux of about 2.3 
        #,iB = 10                         ##<< maximum immobilization flux [1/yr], from Perveen14 Table 1(i), is larger than required immobilization flux of about 2.3 
        ,iB = 25                         ##<< maximum immobilization flux [1/yr], larger than Perveen14 Table 1(i), is larger than required immobilization flux of about 2.3 
        ,iI = 22.91   #0.0627704*365     ##<< input of mineral N [gN/m2/yr] Perveen14 Table 1
        #,l = 2.0075  #0.0055*365    ##<< leaching rate of mineralN l I [1/yr] Perveen14 Table 4
        ,l = 0.95866 #0.00262647 * 365    ##<< leaching rate per mineralN l I [1/gN/yr] Perveen14 Table 1
        #,kIP = 2.19   ##<< plant uptake kIP I [1/yr], here equals plant export eP*CP/betaPlant 0.42*0.0142857*365          
        ,kIP = 16.0351   ##<< plant uptake kIP I [1/yr], here equals plant N input (mp*Cp/CN) + plant export eP*CP/betaPlant 0.00505757*365*525/70  + 0.42*0.0142857*365
        ,useAlphaCUntilNLimited = TRUE      ##<< do not decrease investment into C enzmyes when NSubstrateLimtited, but only when N-Limited
)
parms0 <- within(parms0,{
            kmR <- kmL <- km
            epsR <- epsL <- eps
            cnER <- cnEL <- cnE 
            kNR <- kNL <- kN
        })
parms <- parms0
testScen <- ""

.tmp.testMicrobialLoop <- function(){
	testScen <- "_noTvrMin"
	parms$epsTvr <- parms0$epsTvr <- 1	# nothing of turnover mineralized
	parms$tau <- parms0$tau <- (0.016906/parms0$epsTvr)*365	# also adjust mircobial tvr when changing epsTvr
}

.tmp.f <- function(){
    plantExpN <- 2.19 
    gainCR <- 208.05
    gainNR <- gainCR/parms$cnIR
    #            
    IGain <- with(parms, iI - l - plantExpN)  # 18.71 gN/m2/yr
    (CGain <- IGain*parms$cnIR)
}

# set kR to match A: A*B = kR*R*limER
# parms0$kR <- A_Perveen14 *B00 / (R0 * 0.75)


x0 <- x0Orig <- c( #aE = 0.001*365
        B = B00                 ##<< microbial biomass 
        ,ER  = 2*parms0$km      ##<< total enzyme pool
        ,EL  = 4*parms0$km      ##<< total enzyme pool
        ,R = R0                 ##<< N rich substrate
        ,RN = R0/parms0$cnIR    ##<< N rich substrate N pool
        ,L = L0                 ##<< N poor substrate
        ,LN = L0/parms0$cnIL    ##<< N poor substrate N pool
        ,I =  2.09 -2/5         ##<< inorganic pool g N /m2 Perveeen14, Table 4, start out a little bit lower
)
x <- x0


#pEstNames <- c("kL","iB")
pEstNames <- c("kL","eps")
pEst0 <- pEst <- unlist(parms0[pEstNames])

upperBoundProb = 0.99	# quantile of the upper boundary
parmsBounds = list(		# mode and upper bound
        kL = c(0.5, 10)		
        ,kR = c(1/40, 1/2)
        ,aE = c(0.3, 5)
        #,eps = c(0.4, 0.7)
        ,eps = c(0.3, 0.7)
        ,iB = c(0.1,20)
        ,tau = c(2,20)
)
varDistr <- twVarDistrVec( names(parmsBounds) )	# by default assumed normal
varDistr[c("kL","kR","aE","iB","tau")] <- "lognorm"
varDistr[c("eps")] <- "logitnorm"
parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb, useMedian=FALSE )
poptDistr <- twConstrainPoptDistr(pEstNames, parDistr)
#ggplotDensity.poptDistr(parDistr)

# transform entire parameter vectors between scales
pNorm <- transNormPopt( pEst, parDistr=parDistr )
#transOrigPopt(pNorm,  parDistr=parDistr)

#------------------------ observations

# Perveen14 Table 4
obsOrig <- obs <- c(
     Cf = 635               # fresh litter pool [gC/m2]
    ,N = 2.09               # mineral N pool
    ,leach = 2.0075 #0.0055*365     # N leaching gN/m2/yr
    ,dS = 208.05  #0.57*365          # increase of SOM gC/m2/yr
    ,dR = 208.05             # increase in R pool, for steady state simulation
    #,A = A_Perveen14        # decomposer consumption rate decR/B
    ,immo = 4.871            # immobilization flux to balance inorganic N (not a real observation but assumption that dI=0)
    ,dI = 0                  # assumption of zero change of inorganic N
    ,respO = 0               # assume substrate N limited -> hence no overflow respiration
    ,dB = 0                  # do not match state by non-steady biomass
)
sdObsOrig <- sdObs <- c(
        Cf = 129                # fresh litter pool [gC/m2]
        ,N = 0.68               # mineral N pool
        ,leach = 0.0047*365     # N leaching gN/m2/yr
        ,dS = 0.4*365           # increase of SOM gC/m2/yr
        ,dR = 0.4*365
        #,A = A_Perveen14        # decomposer consumption rate decR/B
        ,immo = as.numeric(obsOrig["immo"]*0.5)     #
        ,dI = as.numeric(obs["N"]/5)           # 0.4g/yr change
        ,respO = 5
        ,dB = 5                 # low change allowed
)

# N fluxes gN/m2/yr
(somAcc <- obs["dR"] / parms$cnB)     #18.91 # N accumulation rate in SOM
(iLN <- parms$iL  / parms$cnIL)    #13.84 # N supply in litter
(leach <- obs["leach"])              #2.0075 leaching
(plantExport <- 7.988e-4*365 *525/parms$cnIL)   # 2.187 plant export eP*CP/cnLitter
(plantUp <- iLN + plantExport)          # 16.031 plant uptake of N
(iI <- parms$iI)                        # 22.910 inorganic N inputs
(immoEq <- iI - leach - plantUp)        # 4.871 immobilization for dI=0
(iLN + immoEq)                          # 18.71 N inputs to SOM nearly equals somAcc 

#------------------ cost function 
costSeam2 <- function(p=parms0, obs=obsOrig, sdObs=sdObsOrig, parDistr=parDistr, parms=parms0, isRecover=FALSE
    #,times = seq(0,100, length.out=301)
    ,times = c(0:6)
    ,x0l=x0
    ,yearsAccR=0       ##<< number of years to decrease R in x0
    ,stopOnError=FALSE  ##<< if lsoda returns error, by default return Inf, set to TRUE to stop
    ,iCosts= c("Cf","dS","respO","dB","dI")
){
    scen <- "Revenue"
    #scen <- "Fixed"
    #resAll <- lapply( c("Revenue","Fixed","Match"), function(scen){
                parms[names(p)] <- pOrig <- transOrigPopt(p, parDistr=parDistr)
                #dur <- 20 
                #times <- c(0,dur,dur+1)
                #times <- seq(0,100, length.out=301)
                #times <- seq(0,10000, length.out=101)
                cnR <- as.numeric(x0l["R"]/x0l["RN"])
                x0Past <- x0l
                x0Past["R"] <- x0l["R"] - yearsAccR*obs["dS"]
                x0Past["RN"] <- x0Past["R"]/cnR
                res <- res1 <- try( as.data.frame(lsoda( x0Past, times, derivSeam2, parms=parms)) )
                #res <- res1 <- as.data.frame(lsoda( x0, times, derivSeam1, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivSeam1, parms=within(parms0, useFixedAlloc<-TRUE) ))
                if( inherits(res,"try-error") ){
                    #return( rep(Inf,5) )
                    if( stopOnError) res <- lsoda( x0Past, times, derivSeam2, parms=parms)
                    return( Inf )
                }
                iCloseI <- which.min( (res[,"I"] - obs["N"])^2 )
                if( isTRUE(parms$isFixedI) ) iCloseI <- nrow(res)
                xE <- unlist(res[iCloseI,])
                dS <- as.numeric(xE["dL"] + xE["dR"])
                #xE[c("PhiB","PhiU","PhiTvr","PhiBU","PhiTotal")]
                pred <-  c(Cf=as.numeric(xE["L"]),N=as.numeric(xE["I"]),leach=as.numeric(xE["I"])*parms$l,dS=dS, dR = as.numeric(xE["dR"])
                , immo=as.numeric(-xE["PhiTotal"]), dI=as.numeric(xE["I"]), respO=as.numeric(xE["respO"]), dB=as.numeric(xE["dB"])
                )
                #relDiffs <- ((obs - pred)/sdObs)
                relDiffs <- ((obs - pred)/sdObs)[iCosts]      # omit N and leaching
                if( isTRUE(isRecover) ){ print("recover in costSeam2"); recover() }
            #})
    ans <- sum(relDiffs^2)
    attr(ans,"relDiffs") <- relDiffs
    attr(ans,"pred") <- pred
    attr(ans,"iCloseI") <- iCloseI
    attr(ans,"resLsoda") <- res
    ans
}

parms00 <- parms0
pNorm0 <- pEst0

.tmp.f <- function(){
    #derivSeam2( 0, x0, within(parms0, isRecover <- TRUE))
    pNorm <- pNorm0
    #pNorm["kL"] <- log(1.6)
    #pNorm["iB"] <- log(2.4)
    pNorm["kL"] <- log(2.6)
    pNorm["eps"] <- logit(0.4)
    pNorm["kR"] <- log(0.054)
    #parms0$isFixedR <- parms0$isFixedI <- TRUE
    #parms0$kR <- 1/70
    x0["B"] <- 140
    parms0$useAlphaCUntilNLimited <- TRUE
    #parms0$useAlphaCUntilNLimited <- FALSE
    (tmp <- costSeam2( pNorm, obsOrig, parDistr=parDistr
                        ,stopOnError=TRUE
                #,isRecover=TRUE
                #,parms=within(parms0, isRecover <- TRUE)
                ))
    as.numeric(tmp <- costSeam2( pNorm, obsOrig, parDistr=parDistr
        , times = seq(0,50, length.out=301)
        #, times = seq(0,25, length.out=301)
        #, times = 1:70
        ,stopOnError=TRUE
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
    ))
    res <- attr(tmp,"resLsoda")
    plotResSeam2(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10"), ylim=c(-100,300))
    #plotResSeam1(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100") )
    #plotResSeam1(res, "topright", cls = c("B10","Rr","I100","MmB100","cnR10"))
    #plotResSeam1(res, "topright", cls = c("cnR"))
    (xE <- as.vector(tail(res,1)[,2:9]))
    #xE <- as.vector(res[56,2:9])
    parmsA <- parms0
    parmsA[names(pNorm)] <- pOrig <- transOrigPopt(pNorm, parDistr=parDistr)
    derivSeam2( 0, xE, within(parmsA, isRecover <- TRUE) )
}

#------------- optimization with constant I and R,L
pEstNames <- c("kL","eps","kR")
pEst0 <- pEst <- unlist(parms0[pEstNames])
poptDistr <- twConstrainPoptDistr(pEstNames, parDistr)
pNorm <- transNormPopt( pEst, parDistr=parDistr )

#parms0$epsTvr <- 0.8
#parms0$tau <- (0.016906/parms0$epsTvr)*365
parmsFixed <- parms0
parmsFixed$isFixedR <- parmsFixed$isFixedI <- parmsFixed$isFixedL <- TRUE
x0Fixed <- x0
x0Fixed["B"] <- 66
iCosts=c("Cf", "immo", "respO", "dB", "dR")
times <- 1:6
tmp <- costSeam2( pNorm, obsOrig, parDistr=parDistr, parms=parmsFixed, x0l=x0Fixed, yearsAcc=0, iCosts=iCosts)
#tmp <- costSeam2( pNorm, obsOrig, parDistr=parDistr, parms=parmsFixed, x0l=x0Fixed, yearsAcc=0, iCosts=iCosts, isRecover=TRUE)
resoFixed <- optim( pNorm, costSeam2, parDistr=parDistr, parms=parmsFixed, control=list(maxit=100L), times=times, iCosts= iCosts, yearsAcc=0 )
pOptFixed <- resoFixed$par
parOptFixed <- transOrigPopt(pOptFixed, parDistr=parDistr)
(tmpOptFixed <- costSeam2( pOptFixed, obsOrig, parDistr=parDistr, parms=parmsFixed, times=times, iCosts=iCosts, yearsAcc=0
                , x0=x0Fixed
                #, x0l=xEOptFixed
                #,isRecover=TRUE
            ))
(xEOptFixed <- (resOptTail <- unlist(tail(attr(tmpOptFixed,"resLsoda"),1)))[2:9])
tmp <- costSeam2( pOptFixed, obsOrig, parDistr=parDistr, parms=parmsFixed, x0l=xEOptFixed, times=seq(0,6, length.out=301), yearsAcc=0)
res <- resOpt <- attr(tmp,"resLsoda")
plotResSeam2(res, "topright", cls = c("B","respO","dR","alpha100","dI","mPhiTotal10"), ylim=c(0,300))

parmsA <- parms0
parmsA[names(pNorm)] <- pOrig <- transOrigPopt(pOptFixed, parDistr=parDistr)
#derivSeam2( 0, xEOptFixed, within(parmsA, isRecover <- TRUE) )

#------------- optimization with changing pools in time
x0LowI <- xEOptFixed;
if( testScen == "_noTvrMin") x0LowI["I"] <- 0.4	# start from low inorganic carbon pools  


iCosts <- c("Cf","dR","respO","dB","N","leach")
as.numeric(tmp <- costSeam2( pOptFixed, obsOrig, parDistr=parDistr, times = seq(0,10, length.out=301)
				, iCosts= iCosts
				, yearsAccR=6, x0=x0LowI
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
        ))
attributes(tmp)$relDiffs
res <-  attr(tmp,"resLsoda")
plotResSeam2(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(0,300))
abline(v=6, col="grey")

# note default iCosts with dS and dI
iCosts <- c("Cf","dR","respO","dB","N","leach")
reso <- optim( pOptFixed, costSeam2, parDistr=parDistr, control=list(maxit=100)
        , iCosts= iCosts
        , yearsAccR=6, x0=x0LowI
)
pOpt <- reso$par
parOpt <- transOrigPopt(pOpt, parDistr=parDistr)
(tmpOpt <- costSeam2( pOpt, obsOrig, parDistr=parDistr, yearsAccR=6, x0=x0LowI, iCosts=iCosts))
(xEOpt <- (resTailOpt <- unlist(tail(attr(tmpOpt,"resLsoda"),1))[2:9]))

as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = seq(0,10, length.out=301)
        ,yearsAccR=6, x0=x0LowI
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
        ))
iClose <-  attr(tmp,"iCloseI")
#if( testScen == "_noTvrMin") iClose <- 30	# before L and R are not in steady state  
res <- resOpt50 <- attr(tmp,"resLsoda")
plotResSeam2(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(0,300))
timeClose <- resOpt50$time[iClose]
abline(v=timeClose, col="grey")

#derivSeam2( 0, xEOpt, within(parmsA, isRecover <- TRUE))

# take state after fluctuations as Base for experiments
xBaseSt <- unlist(resOpt50[iClose,2:9])
timesExp = seq(0,50, length.out=501)
timesExpFit <- 0

as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = timesExp
                ,x0l=xBaseSt
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
        ))
res <- resControl <- attr(tmp,"resLsoda")
plotResSeam2(res[res$time<=10,], "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS","mPhiBU10"), ylim=c(0,300))
abline(v=timesExpFit, col="grey")

(xE <- unlist(res[80,2:9]))
#xE <- as.vector(res[56,2:9])
parmsA <- parms0
parmsA[names(pNorm)] <- pOrig <- transOrigPopt(pOpt, parDistr=parDistr)
#derivSeam2( 0, xE, within(parmsA, isRecover <- TRUE) )

# when simulating a longer period, N-substrate limited system
# inorganic pool increases (need to include plant model to change this behaviour)
# R accumulation slowly declines

.tmp.f <- function(){
    x0IncI <- xBaseSt
    x0IncI["I"] <- 1
    parmsCLim <- within(parms0,{
                useAlphaCUntilNLimited <- TRUE
            })
    as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = timesExp
                    ,parms=parmsCLim
                    ,x0l=x0IncI
                    ,yearsAccR=-2
            #,isRecover=TRUE
            #,parms=within(parms0, isRecover <- TRUE)
            ))
    res <- attr(tmp,"resLsoda")
    plotResSeam2(res[res$time<=10,], "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(0,300))
    abline(v=timesExpFit, col="grey")
    (xETmp <- unlist(tail(res,))[2:9])
    (xETmp <- unlist(res[which.min((res$time-8)^2),][2:9]))
    parmsA <- parms0
    parmsA[names(pNorm)] <- pOrig <- transOrigPopt(pOpt, parDistr=parDistr)
    #derivSeam2( 0, xETmp, within(parmsA, isRecover <- TRUE) )
    
    #cnR <- xETmp["R"]/xETmp["RN"] 
    #xETmp["R"] <- 12000
    #xETmp["RN"] <- xETmp["R"] / cnR
    as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = timesExp
                    ,parms=parmsCLim
                    ,x0l=xETmp
            #,isRecover=TRUE
            #,parms=within(parms0, isRecover <- TRUE)
            ))
    res <-  attr(tmp,"resLsoda")
    plotResSeam2(res[res$time<=10,], "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(0,300))
    
    xBaseSt <- xETmp
    parms0 <- parmsCLim
}

#----------- experiment: increased plant C input 
# C input by plants increased by 6/4 and increased C/N input from 70 to 90
parmsEl <- within(parms0,{
            iL <- parms0$iL*6/4
            cnIL <- parms0$cnIL*5/4
            # tau <- tau / 1.2 does not change R accumulation nor enzyme partitioning (alpha)
        })
(iLN <- parmsEl$iL  / parmsEl$cnIL)    #13.84 # N supply in litter
(plantExport <- 7.988e-4*365 *525/parms$cnIL)   # 2.187 plant export eP*CP/cnLitter
(plantUp <- iLN + plantExport)          # 16.031 plant uptake of N
parmsEl$kIP <- plantUp


as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = timesExp
                ,x0l=xBaseSt
                ,stopOnError=TRUE
                ,parms=parmsEl
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
        ))
resEl <- res <- attr(tmp,"resLsoda")
plotResSeam2(resEl, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(0,300))
plotResSeam2(resEl[resEl$time<=10,], "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(0,300))

# becomes N limited -> C overflow
# shift towards degrading R pool
# fast increase of L to a new level
# together with increased microbial biomass and biomass turnover, still accumulation in R at the same rate
# uses N inputs to build SOM, decreased leaching

#-------------- experiment: decreased N inputs
parmsIlow <- within(parms0,{
            iI <- 1  #0.01*365      # still larger than leaching losses -> N Accumulation
            #iI <- 0  #0.01*365
            cnIL <- cnIL*2
        })
# assuming the plant still takes up the what supplied by litter
(iLN <- parmsIlow$iL  / parmsIlow$cnIL)    #13.84 # N supply in litter
(plantExport <- 7.988e-4*365 *525/parms$cnIL)   # 2.187 plant export eP*CP/cnLitter
(plantUp <- iLN + plantExport)          # 16.031 plant uptake of N
parmsIlow$kIP <- plantUp


as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = timesExp
                ,yearsAccR=0, x0l=xBaseSt
                ,stopOnError=TRUE
                ,parms=parmsIlow
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
        ))
resILow <- res <-attr(tmp,"resLsoda")
#plotResSeam1(resIlow, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(0,300))
plotResSeam2(resILow, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(-100,300))
plotResSeam2(resILow[resILow$time<=10,], "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(-100,300))
abline(h=0, col="grey")


tail(res)

(xE <- as.vector(tail(res,1)[,2:9]))
#xE <- as.vector(res[56,2:9])
parmsA <- parms0
parmsA[names(pNorm)] <- pOrig <- transOrigPopt(pOpt, parDistr=parDistr)
#derivSeam2( 0, xE, within(parmsA, isRecover <- TRUE) )

# Inorganic pool gets so low, that microbes cannot compensate their uptake efficiency losses by 
# immobilization -> availability for plants

#-------------- experiment: increased N inputs
parmsIHigh <- within(parms0,{
            iI <- 25.6  #0.07*365      # still larger than leaching losses -> N Accumulation
            #cnIL <- cnIL*0.9
        })
# assuming the plant still takes up the what supplied by litter
(iLN <- parmsIHigh$iL  / parmsIHigh$cnIL)    #13.84 # N supply in litter
(plantExport <- 7.988e-4*365 *525/parms$cnIL)   # 2.187 plant export eP*CP/cnLitter
(plantUp <- iLN + plantExport)          # 16.031 plant uptake of N
parmsIHigh$kIP <- plantUp


as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = timesExp
                ,yearsAccR=0, x0l=xBaseSt
                ,stopOnError=TRUE
                ,parms=parmsIHigh
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
        ))
resIHigh <- res <-attr(tmp,"resLsoda")
plotResSeam2(resIHigh, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(-100,300))
plotResSeam2(resIHigh[resIHigh$time < 15,], "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","mPhiTotal10","dR","dS"), ylim=c(-100,300))


tail(res)

(xE <- unlist(res[80,2:9]))
#xE <- as.vector(res[56,2:9])
parmsA <- parms0
parmsA[names(pNorm)] <- pOrig <- transOrigPopt(pOpt, parDistr=parDistr)
#derivSeam2( 0, xE, within(parmsA, isRecover <- TRUE) )

#---------------- plotting
resScens <- rbind(
        cbind( scenario="control", resControl)
        ,cbind( scenario="elevated CO2", resEl)
        ,cbind( scenario="decreased N inputs", resILow)
        ,cbind( scenario="increased N inputs", resIHigh)
)
resScens$dS <- resScens$dR + resScens$dL
resScens$leaching <- resScens$I * parms0$l
predM <- melt(resScens, 1:2)
predM$variable <- relevel(relevel(relevel( relevel(predM$variable, "dR"),"R"),"L"),"alpha")
levels(predM$variable)[match(c("L","R","dR","I","leaching","PhiBU","PhiTotal"), levels(predM$variable))] <- 
		c("L~(gC~m^{-2})","R~(gC~m^{-2})","dR~(gC~m^{-2}*yr^{-1})","I~(gN~m^{-2})","leaching~(gN~m^{-2}*yr^{-1})","Phi_BU~(gN~m^{-2}*yr^{-1})","Phi~(gN~m^{-2}*yr^{-1})")
predMCtrl <- subset(predM, scenario %in% "control")

dsObs <- data.frame(value=obs, sd=sdObs, time=timesExpFit, scenario="control")
dsObs$variable <- relevel(relevel( relevel(as.factor(rownames(dsObs)), "N"), "dR"),"Cf")
levels(dsObs$variable)[match(c("Cf","dR","N","leach"), levels(dsObs$variable))] <- c("L~(gC~m^{-2})","dR~(gC~m^{-2}*yr^{-1})","I~(gN~m^{-2})","leaching~(gN~m^{-2}*yr^{-1})")
#dssObs <- subset( dsObs, dsObs$variable %in% dsObs$variable[na.omit(pmatch(c("alpha","L","dR","I","leaching"), dsObs$variable))] )
dssObs <- dsObs[grep(c("^alpha|^L|^dR|^I|^leach"),dsObs$variable),]
#p1 <- ggplot( dss <- subset( predMCtrl[grep(c("^alpha$|^L~|^dR~|^I~|^Phi~|^Phi~|^leach"),predMCtrl$variable),], time <=5)
#p1 <- ggplot( dss <- subset( predMCtrl[grep(c("^alpha$|^L~|^dR~|^I~|^Phi~|^leach"),predMCtrl$variable),], time <=5)
p1 <- ggplot( dss <- subset( predMCtrl[grep(c("^L~|^dR~|^I~|^leach"),predMCtrl$variable),], time <=5)
        , aes(x=time, y=value, linetype=scenario, colour=scenario)) +
        geom_line(size=baseLineSize) +
        geom_point( data=dssObs, colour="black" )+
        geom_errorbar( data=dssObs,aes(ymin=value-sd, ymax=value+sd), width=.25, colour="black") +
        #facet_grid( variable ~ .,  scales="free_y", labeller = label_parsed) +
        #facet_wrap( ~ variable,  scales="free_y", labeller = label_parsed) +
		facet_wrap( ~ variable,  scales="free_y", labeller = label_parsed) +
		themeDefault + 
		theme(strip.background = element_blank()) +			
        #theme_bw(base_size=baseFontSize) +
        theme(axis.title.y = element_blank()) + 
        xlab("time (yr)") +
        theme(legend.position = "none") +
		annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
        c()
twWin(3.27, 2)
p1
savePlot(paste0("soilPaper14/fig/pastureFitMatch",testScen,".pdf"),"pdf")

#library(ggthemes) # theme_tufte does not display axes lines
p1b <- ggplot( dss <- subset( predM[grep(c("^alpha$|^L~|^dR~|^I~|^Phi~"),predM$variable),], time <=5)
        , aes(x=time, y=value, linetype=scenario, colour=scenario)) +
		geom_hline(yintercept=0, color="grey75") +		
		geom_line(size=baseLineSize) +
        #facet_wrap( ~ variable,  scales="free_y", nrow=6, ncol=1) +
        #facet_wrap( ~ variable,  scales="free_y", ncol=2L, labeller = label_parsed) +
		facet_wrap( ~ variable,  scales="free", ncol=2L, labeller = label_parsed) +
		themeDefault +
		#theme_tufte() +
		theme(strip.background = element_blank()) +			
		theme(axis.title.y = element_blank()) + 
        xlab("time (yr)") +
		#annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
		#theme(legend.position = "bottom") +
        theme(legend.position = c(0.95,-0.08), legend.justification=c(1,0)) +
        theme(legend.title = element_blank()) +
        guides(linetype=guide_legend(nrow=4,byrow=TRUE)) +
		#theme(panel.grid.major.x=element_line(colour="grey75")) +
		#annotate("segment", x=1:5, xend=1:5, y=-Inf, yend=)
        c()
twWin(3.27, 3.5)
p1b
# use scale="free" to get x-axis ticks
# then remove lables
#http://stackoverflow.com/questions/38986702/add-tick-marks-to-facet-plots-in-r
g1<-ggplotGrob(p1b)
#g1
# look for axis-b and dig down what are the lables
for( g in c(16,17,18,19) )
	g1[[1]][[g]]$children[[2]]$grobs[[2]] <- NULL #nullGrob()
grid.newpage()
grid.draw(g1)
#savePlot("soilPaper14/fig/pastureFitInputScenarios.pdf","pdf")


twWin(3.27, 5.9)
#twWin(7.5,8.7)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p1b, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
#savePlot(paste0("soilPaper14/fig/pastureFit",testScen,".pdf","pdf")



#------------- plotting differently for presentations
.tmp.f <- function(){
    twWin(7,4.6)
    p1p <- ggplot( dss <- subset( predMCtrl[grep(c("^L~|^R~|^dR~|^I~|^Phi~|^leach"),predMCtrl$variable),], time <=5)
                    , aes(x=time, y=value, linetype=scenario, colour=scenario)) +
            geom_line(size=baseLineSize) +
            geom_point( data=dssObs, colour="black", size=3 )+
            geom_errorbar( data=dssObs,aes(ymin=value-sd, ymax=value+sd), width=.25, colour="black", size=baseLineSize) +
            facet_wrap( ~ variable,  scales="free_y", labeller = label_parsed) +
            theme_bw(base_size=baseFontSize) +
            theme(axis.title.y = element_blank()) + 
            xlab("time (yr)") +
            theme(legend.position = "none") +
            theme()
    p1p
    
    twWin(5.7,4.6)
    p1bp <- ggplot( dss <- subset( predM[grep(c("^alpha$|^L~|^dR~|^Phi~|^I~"),predM$variable),], time <=5)
                    , aes(x=time, y=value, linetype=scenario, colour=scenario)) +
            geom_line(size=baseLineSize) +
            facet_wrap( ~ variable,  scales="free_y", labeller = label_parsed) +
            #facet_grid( variable ~ .,  scales="free_y") +
            theme_bw(base_size=baseFontSize) +
            theme(axis.title.y = element_blank()) + 
            xlab("time (yr)") +
            theme(legend.position = "bottom") +
            theme(legend.title = element_blank()) +
            guides(linetype=guide_legend(nrow=2,byrow=TRUE)) +
            theme()
    p1bp
    direct.label(p1bp, list(last.points, hjust = 0.7, vjust = 1))
}    

# intitial state for table
xBaseSt