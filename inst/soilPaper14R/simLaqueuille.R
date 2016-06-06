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
        ,kN= 1/(1/12)   ##<< enzyme tvr 1 month (Blagodatskaya 60 days priming)
        ,kNB = 0.8      ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R), fixed
        ,kR = 1/(10)    ##<< 1/(x years) 
        #,kR = A_Perveen14 * B00 / (R0 *0.5)  #1/(10)     ##<< 1/(x years)     estimated from A in Perveen14
        #,kL = 2.4 #1/(0.33)   ##<< 1/(x years)     # formerly 1 year  # 4.23e-4*365
        ,kL = 2 # 1/yr, turns over X per year 
        #,aE = 1e-4*365 #0.001*365  ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,aE = 1e-3*365   ##<< C biomass allocated to enzymes gC/yr/microbial biomass, corresponds about 2% of uptake, toCheck
        ,km = 0.08       ##<< enzyme half-saturation constant
        #,km = 0.03     ##<< enzyme half-saturation constant
        ,m = 0 # 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = (0.016906)*365 #1/60*365  ##<< biomass turnover rate from Perveen (s)
        #,tau = 1/10*365 #      ##<< biomass turnover rate of x days, with too low turnover time, too low microbial biomass
        #,tau = (0.016906+0.0368857)*365 #1/60*365  ##<< biomass turnover rate from Perveen (s+r) to match similar microbes pool
        #,tau = 1/20*365 #      ##<< biomass turnover rate of x days 
        ,eps = 0.33        ##<< carbon use efficiency  Perveen14: r/(s+r) 31% (here need to addionally account for 2% for enzymes)
        ,epsTvr = 1 #0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
        ,iR = 0         ##<< input by feedback modelled explicitely
        ,iL = 0.00505757*365*525         # g/m2 input per year mp*Cp
        ,plantNUp = 0                    # rate of organic N uptake
        ,useFixedAlloc=FALSE    ##<< set to true to use Fixed enzyme allocation (alpha = 0.5)
        #
        #,kIP = 10.57 #0.0289652*365          ##<< plant uptake kIP I [1/yr]        
        ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I [1/yr]
        ,iI = 22.91   #0.0627704*365     ##<< input of mineral N [gN/m2/yr] Perveen14 Table 1
        #,l = 2.0075  #0.0055*365    ##<< leaching rate of mineralN l I [1/yr] Perveen14 Table 4
        ,l = 0.95866 #0.00262647 * 365    ##<< leaching rate per mineralN l I [1/gN/yr] Perveen14 Table 1
        #,kIP = 2.19   ##<< plant uptake kIP I [1/yr], here equals plant export eP*CP/betaPlant 0.42*0.0142857*365          
        ,kIP = 16.0351   ##<< plant uptake kIP I [1/yr], here equals plant N input (mp*Cp/CN) + plant export eP*CP/betaPlant 0.00505757*365*525/70  + 0.42*0.0142857*365          
)
parms0 <- within(parms0,{
            kmR <- kmL <- km
            epsR <- epsL <- eps
            cnER <- cnEL <- cnE 
            kNR <- kNL <- kN
        })
parms <- parms0

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
        ,I =  2.09              ##<< inorganic pool g N /m2 Perveeen14, Table 4
)
x <- x0


pEstNames <- c("kL","iB")
pEst0 <- pEst <- unlist(parms0[pEstNames])

upperBoundProb = 0.99	# quantile of the upper boundary
parmsBounds = list(		# mode and upper bound
        kL = c(0.5, 10)		
        ,kR = c(1/40, 1/2)
        ,aE = c(0.3, 5)
        #,eps = c(0.4, 0.7)
        ,eps = c(0.3, 0.7)
        ,iB = c(0.1,20)
)
varDistr <- twVarDistrVec( names(parmsBounds) )	# by default assumed normal
varDistr[c("kL","kR","aE","iB")] <- "lognorm"
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
    ,dR = 208.05  #0.57*365          # increase of SOM gC/m2/yr
    #,A = A_Perveen14        # decomposer consumption rate decR/B
)
sdObsOrig <- sdObs <- c(
        Cf = 129                # fresh litter pool [gC/m2]
        ,N = 0.68               # mineral N pool
        ,leach = 0.0047*365     # N leaching gN/m2/yr
        ,dR = 0.4*365           # increase of SOM gC/m2/yr
        #,A = A_Perveen14        # decomposer consumption rate decR/B
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
    ,times = c(0,20,21)
    ,x0l=x0
){
    scen <- "Revenue"
    #scen <- "Fixed"
    #resAll <- lapply( c("Revenue","Fixed","Match"), function(scen){
                parms[names(p)] <- pOrig <- transOrigPopt(p, parDistr=parDistr)
                #dur <- 20 
                #times <- c(0,dur,dur+1)
                #times <- seq(0,100, length.out=301)
                #times <- seq(0,10000, length.out=101)
                cnR <- x0l["R"]/x0l["RN"]
                x0Past <- x0l
                x0Past["R"] <- x0l["R"] - 20*obs["dR"]
                x0Past["RN"] <- x0Past["R"]/cnR
                res <- res1 <- try( as.data.frame(lsoda( x0Past, times, derivSeam2, parms=parms)) )
                #res <- res1 <- as.data.frame(lsoda( x0, times, derivSeam1, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivSeam1, parms=within(parms0, useFixedAlloc<-TRUE) ))
                if( inherits(res,"try-error") ){
                    #return( rep(Inf,5) )
                    return( Inf )
                }
                xE <- unlist(tail(res,1))
                xE2 <- tail(res,2)
                dR <- diff(xE2[,"R"]) / diff(xE2[,"time"])
                ASim <- parms$kR * R0 * xE["limER"] / xE["B"]
                #pred <-  c(Cf=xE["L"],N=xE["I"],leach=xE["I"]*parms$l,dR=dR,A=ASim)
                pred <-  c(Cf=as.numeric(xE["L"]),N=as.numeric(xE["I"]),leach=as.numeric(xE["I"])*parms$l,dR=dR)
                relDiffs <- ((obs - pred)/sdObs)
                if( isTRUE(isRecover) ){ print("recover in costSeam2"); recover() }
            #})
    ans <- sum(relDiffs^2)
    attr(ans,"relDiffs") <- relDiffs
    attr(ans,"pred") <- pred
    attr(ans,"resLsoda") <- res
    ans
}

parms00 <- parms0
pNorm0 <- pEst0

.tmp.f <- function(){
    #derivSeam2( 0, x0, within(parms0, isRecover <- TRUE))
    pNorm <- pNorm0
    pNorm["kL"] <- log(1.6)
    pNorm["iB"] <- log(2.4)
    parms0$kR <- 1/70
    x0["B"] <- 50
    as.numeric(tmp <- costSeam2( pNorm, obsOrig, parDistr=parDistr
        , times = seq(0,25, length.out=301)
        #, times = 1:70
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
    ))
    res <- attr(tmp,"resLsoda")
    plotResSeam1(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100"), ylim=c(0,300))
    #plotResSeam1(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100") )
    #plotResSeam1(res, "topright", cls = c("B10","Rr","I100","MmB100","cnR10"))
    #plotResSeam1(res, "topright", cls = c("cnR"))
    (xE <- as.vector(tail(res,1)[,2:9]))
    #xE <- as.vector(res[56,2:9])
    parmsA <- parms0
    parmsA[names(pNorm)] <- pOrig <- transOrigPopt(pNorm, parDistr=parDistr)
    derivSeam2( 0, xE, within(parmsA, isRecover <- TRUE) )
}


#------------- optimization
(tmp <- costSeam2( pNorm, obsOrig, parDistr=parDistr))
reso <- optim( pNorm, costSeam2, parDistr=parDistr)
pOpt <- reso$par
(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr))
attr(tmp,"relDiffs")
res <- resOpt <- attr(tmp,"resLsoda")
#tail(resOpt,1)
(xEOpt <- as.vector(tail(attr(tmp,"resLsoda"),1)[,2:9]))

as.numeric(tmp <- costSeam2( pOpt, obsOrig, parDistr=parDistr, times = seq(0,50, length.out=301)
        #,isRecover=TRUE
        #,parms=within(parms0, isRecover <- TRUE)
        ))
res <- resOpt50 <- attr(tmp,"resLsoda")
plotResSeam1(resOpt50, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","LN"), ylim=c(0,300))

# when simulating a longer period, N-substrate limited system becomes C-limited
# then inorganic pool increases (need to include plant model to change this behaviour)
# R accumulation does not change

#----------- experiments 
# C input by plants increased by 6/4 and increased C/N input from 70 to 90
parmsEl <- within(parms0,{
            iL <- iL*6/4
            cnIL <- cnIL*7/4
            # tau <- tau / 1.2 does not change R accumulation nor enzyme partitioning (alpha)
        })
(iLN <- parmsEl$iL  / parmsEl$cnIL)    #13.84 # N supply in litter
(plantExport <- 7.988e-4*365 *525/parms$cnIL)   # 2.187 plant export eP*CP/cnLitter
(plantUp <- iLN + plantExport)          # 16.031 plant uptake of N
#parmsEl$kIP <- plantUp

x0El <- x0
#x0El["ER"] <- 0.52
#x0El["EL"] <- 0.1

resEl <- res <- attr(tmp <- costSeam2( pOpt, obsOrig, parms=parmsEl, parDistr=parDistr, times = seq(0,25, length.out=301)
    ,x0=x0El
),"resLsoda")
plotResSeam1(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100","LN"), ylim=c(0,300))
plotResSeam1(res, "topright", cls = c("B","respO","Rr","Lr","alpha100","I100") )
(xEel <- as.vector(tail(resEl,1)[,2:9]))
#xE <- as.vector(res[56,2:9])
parmsA <- parms0
parmsA[names(pOpt)] <- pOrig <- transOrigPopt(pOpt, parDistr=parDistr)
derivSeam2( 0, xEel, within(parmsA, isRecover <- TRUE) )

# N content in litter does change only slightly, enzyme partitioning does not change -> all to overflow 
