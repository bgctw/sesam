quit()
isPaperBGC <- ("paperBGC" %in% commandArgs(trailingOnly = TRUE))
#isPaperBGC <- TRUE

# simulating CO2 increase and bare soil, based on modEezy5
baseFontSize <- 16  # presentations
# based on 
if( isPaperBGC){
    library(twDev)
    loadPkg()
    baseFontSize <- 9  # pubs
}
library(ggplot2)
library(grid)   #unit
library(reshape)  # melt
library(RColorBrewer) # brewer.pal
library(twDEMC)
library(twDEMCPlot)
library(FME)

# gC/m2 and gN/m2, /yr
# From Perveen 2014 Symphony
B00 <- 16
A_Perveen14 <- 0.0317917*365 
L0 <- L0_Perveen14 <- 635
SOM0_Allard07 <- 12000
R0 <-  SOM0_Allard07 - L0_Perveen14   #640*0.1*0.18*1000 rho*10cm*Corg

NLoss_A_Perveen14 <- 0.0055 # phi_i - l*N - phi_up

parms0 <- list(
        cnB = 11        ##<< Perveen14
        ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        ,cnIR = 4.5    ##<< between micr and enzyme signal (not used)
        ,cnIL = 70      ##<< N poor substrate, Perveen14
        #,kN= 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
        ,kN= 1/(1/12)   ##<< enzyme tvr 1 month (Blagodatskaya 60 days priming)
        ,kNB = 0.8      ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R), fixed
        #,kR = 1/(50)    ##<< 1/(x years) 
        ,kR = A_Perveen14 * B00 / (R0 *0.5)  #1/(10)     ##<< 1/(x years)     estimated from A in Perveen14
        ,kL = 2.4 #1/(0.33)   ##<< 1/(x years)     # formerly 1 year  # 4.23e-4*365
        ,aE = 1e-4*365 #0.001*365  ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,km = 0.08       ##<< enzyme half-saturation constant
        #,km = 0.03     ##<< enzyme half-saturation constant
        ,m = 0 # 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = (0.016906+0.0368857)*365 #1/60*365  ##<< biomass turnover rate
        ,eps = 0.35 #0.3       ##<< carbon use efficiency  Perveen14: 30% (here +tvrB and tvrE)
        ,epsTvr = 1 #0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
        ,iR = 0         ##<< input by feedback modelled explicitely
        ,iL = 0.00505757*365*525         # g/m2 input per year mp*Cp
        ,plantNUp = 0
        ,useFixedAlloc=FALSE    ##<< set to true to use Fixed enzyme allocation (alpha = 0.5)
        #
        ,iP = 10.57 #0.0289652*365          ##<< plant uptake iP I [1/yr]        
        ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I [1/yr]
        ,iI = 23   #0.063*365     ##<< input of mineral N [gN/m2/yr]
        ,l = 0.96   #0.00262647*365       ##<< leaching rate of mineralN l I [1/yr]
)
parms0 <- within(parms0,{
            kmR <- kmL <- km
            epsR <- epsL <- eps
            cnER <- cnEL <- cnE 
            kNR <- kNL <- kN
        })
parms <- parms0

# set kR to match A: A*B = kR*R*limER
parms0$kR <- A_Perveen14 *B00 / (R0 * 0.75)

x0 <- x0Orig <- c( #aE = 0.001*365
        B = 17                     ##<< microbial biomass 
        ,ER  = 2*parms0$km                  ##<< total enzyme pool
        ,EL  = 4*parms0$km                   ##<< total enzyme pool
        ,R = R0                 ##<< N rich substrate
        ,RN = R0/parms0$cnIR    ##<< N rich substrate N pool
        ,L = L0                  ##<< N poor substrate
        ,LN = L0/parms0$cnIL    ##<< N poor substrate N pool
        ,I =  0                    ##<< inorganic pool
)
x <- x0



obsOrig <- obs <- c(
     Cf = 635       # fresh litter pool [gC/m2]
    ,N = 2.09
    ,leach = 0.0055*365     # N leaching gN/m2/yr
    ,dR = 0.57*365   # increase of SOM gC/m2/yr
    ,A = A_Perveen14    # decomposer consumption rate decR/B
)

costSeam2 <- function(p=parms0, obs=obsOrig, parDistr=parDistr, parms=parms0, isRecover=FALSE, useRelativeError=FALSE
){
    scen <- "Revenue"
    #scen <- "Fixed"
    #resAll <- lapply( c("Revenue","Fixed","Match"), function(scen){
                parms[names(p)] <- pOrig <- transOrigPopt(p, parDistr=parDistr)
                dur <- 20 
                times <- c(0,dur,dur+1)
                #times <- seq(0,100, length.out=301)
                #times <- seq(0,10000, length.out=101)
                x0Past <- x0
                x0Past["R"] <- x0["R"] - 20*obs["dR"]            
                res <- res1 <- try( as.data.frame(lsoda( x0Past, times, derivSeam2, parms=parms)) )
                #res <- res1 <- as.data.frame(lsoda( x0, times, derivSeam1, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivSeam1, parms=within(parms0, useFixedAlloc<-TRUE) ))
                if( inherits(res,"try-error") ){
                    return( rep(Inf,5) )
                }
                xE <- unlist(tail(res,1))
                xE2 <- tail(res,2)
                dR <- diff(xE2[,"R"]) / diff(xE2[,"time"])
                ASim <- parms$kR * R0 * xE["limER"] / xE["B"]
                pred <-  c(xE["L"],xE["I"],xE["I"]*parms$l,dR,ASim)
                resids <- -obs + pred
                if( isTRUE(useRelativeError) )
                    resids <- (-obs + pred)/obs
                if( isTRUE(isRecover) ){ print("recover in costSeam2"); recover() }
            #})
    resids
}
.tmp.f <- function(){
    resids <- c(
            L = as.numeric((xE["L"] - obs["Cf"])/obs["Cf"])
            ,N = as.numeric((x["I"] - obs["N"])/obs["N"])
            ,lN = as.numeric((x["I"]*parms$l - obs["leach"])/obs["leach"])
            ,dR = as.numeric((dR - obs["dR"])/obs["dR"])
            ,A = as.numeric((ASim - obs["A"] )/obs["A"])
    )
    derivSeam2( 0, x=xE[1:length(x0Past)+1])
    plotResSeam1(res, "topright", cls = c("B10","respO","Rr","Lr","alpha100"))
}



#pEstNames <- c("kL","kR","aE","eps")
pEstNames <- c("kL","kR","aE")
#pEstNames <- c("kL","kR","eps")
pEst0 <- pEst <- unlist(parms0[pEstNames])

upperBoundProb = 0.99	# quantile of the upper boundary
parmsBounds = list(		# mode and upper bound
        kL = c(0.5, 2)		
        ,kR = c(1/40, 1/2)
        ,aE = c(0.3, 5)
        #,eps = c(0.4, 0.7)
        ,eps = c(0.3, 0.7)
)
varDistr <- twVarDistrVec( names(parmsBounds) )	# by default assumed normal
varDistr[c("kL","kR","aE")] <- "lognorm"
varDistr[c("eps")] <- "logitnorm"
parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb, useMedian=FALSE )
poptDistr <- twConstrainPoptDistr(pEstNames, parDistr)
#ggplotDensity.poptDistr(parDistr)

# transform entire parameter vectors between scales
pNorm <- transNormPopt( pEst, parDistr=parDistr )
#transOrigPopt(pNorm,  parDistr=parDistr)

costSeam2( pNorm, obsOrig, parDistr=parDistr )
#costSeam2( pNorm, obsOrig, parDistr=parDistr, isRecover=TRUE )


#library(FME)
#(resFit <- modFit( costSeam2, pNorm, obs=obsOrig, parDistr=parDistr )) #takes long
#(resFit <- modFit( costSeam2, pNorm, obs=obsOrig, parDistr=parDistr, isRecover=TRUE )) #takes long

parmsi <- parms0

pEstNorm <- resFit$par
costSeam2( pEstNorm, obsOrig, parDistr=parDistr )
(pEst <- parms <- unlist(transOrigPopt(pEstNorm, parDistr=parDistr)))
parms <- parmsi; parms[names(pEst)] <- pEst
times <- seq(0,100, length.out=301)
x0Past <- x0
x0Past["R"] <- x0["R"] - 20*obs["dR"]            
res <- res1 <- as.data.frame(lsoda( x0Past, times, derivSeam2, parms=parms))
xE <- unlist(tail(res,1))
plotResSeam1(res, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
plotResSeam1(res, "topright", cls = c("B10","Mm","Rr","Lr","alpha100"))
abline(v=20, lty="dashed", col="grey")

#pEstNorm1 <- pEstNorm

#----- continue with including eps

pEstNames <- c("kL","kR","aE","eps")
pNorm <- c( pEstNorm1, transNormPopt( unlist(parms0["eps"]), parDistr=parDistr ))
#transOrigPopt(pNorm,  parDistr=parDistr)

costSeam2( pNorm, obsOrig, parDistr=parDistr )
#costSeam2( pNorm, obsOrig, parDistr=parDistr, isRecover=TRUE )


#library(FME)
#(resFit <- modFit( costSeam2, pNorm, obs=obsOrig, parDistr=parDistr )) #takes long

# local optimum for given eps

#---------- prescribe different values of eps and optimize for each
pEstNames <- c("kL","kR","aE")
#pEstNames <- c("kL","kR","eps")
pEst0 <- pEst <- unlist(parms0[pEstNames])
pNorm <- transNormPopt( pEst, parDistr=parDistr )

epsis <- seq(0.3,0.37,by=0.002)
resLeps <- resL <- lapply( epsis, function(epsi){
            cat(epsi,",")
            parmsi <- within(parms0, eps <- epsi)
            resFit <- modFit( costSeam2, pNorm, obs=obsOrig, parDistr=parDistr, parms=parmsi )            
        })
names(resL) <- epsis
ssr <- sapply( resL, "[[", "ssr")
ssr[ssr > 40] <- NA
plot( ssr ~ epsis )

iBest <- which.min(ssr)
resFitAbs <- resFit <- resL[[iBest]]
epsBest <- epsis[iBest]

#------ repeat with relative errors
pNorm <- resFitAbs$par
parmsi <- within(parms0, eps <- epsBest)
resFitRel <- resFit <- modFit( costSeam2, pNorm, obs=obsOrig, parDistr=parDistr, parms=parmsi, useRelativeError=TRUE )            

.tmp.f <- function(){
    #------ repeat with relative errors for standard esp=0.3
    epsi <- 0.3
    pNorm <- resLeps[[1]]$par
    parmsi <- within(parms0, eps <- epsi)
    resFitRel03 <- resFit <- modFit( costSeam2, pNorm, obs=obsOrig, parDistr=parDistr, parms=parmsi, useRelativeError=TRUE )            
    # no decompostion of R
}



