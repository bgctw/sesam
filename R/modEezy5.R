#library(deSolve)
# feeding back microbial turnover and enzymes to the S1 pool, 
# now explicitly modelling N because CN-ratio can change
# 

# gC/m2 and gN/m2, /yr
parms0 <- list(
        cnB = 7.16
        ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        #,cnE = 7.16
        ,cnIS1 = 4.5      ##<< between micr and enzyme signal
        ,cnIS2 = 30 ##<< N poor substrate
        #,kN = 0.05   ##<< (per day) enzyme turnover
        ,kN= 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
        ,kNB = 0.8     ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of S1)
        #,kS1 = 0.2      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        #,kS2 = 1      ##<< substrate decomposition rate N-poor
        #,kS1 = 5e-3      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        #,kS2 = 10e-3     ##<< substrate decomposition rate N-poor
        #,kS2 = 5e-2     ##<< substrate decomposition rate N-poor
        #,aE = 0.05   ##<< C-uptake allocated to enzymes
        #,kS1 = 1/(5*365)        ##<< 5 years 
        #,kS2 = 1/(0.5*365)        ##<< 1/2 years 
        ,kS1 = 1/(50)        ##<< 1/(x years) 
        ,kS2 = 1/(1)        ##<< 1/(x years) 
        #,aE = 0.003*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,aE = 0.001*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,km = 0.3     ##<< enzyme half-saturation constant
        #,km = 0.03     ##<< enzyme half-saturation constant
        #,km = 14     ##<< enzyme half-saturation constant
        ,m = 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = 1/60*365  ##<< biomass turnover rate (12 days)
        ,eps = 0.5      ##<< carbon use efficiency
        ,epsTvr = 0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
        ,iS1 = 0        ##<< input modelled explicitely
        ,iS2 = 300         # g/m2 input per year (half NPP)
        #,plantNUp = 300/70*1/4  # plant N uptake balancing N inputs
        ,plantNUp = 0
        ,useFixedAlloc=FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
)
parms0 <- within(parms0,{
 km1 <- km2 <- km
 eps1 <- eps2 <- eps
 cnE1 <- cnE2 <- cnE 
 kN1 <- kN2 <- kN
        })

parms <- parms0

x0 <- x0Orig <- c(  # aE = 0.003*365
        B = 1                     ##<< microbial biomass 
        ,E1  = 2*parms0$km                  ##<< total enzyme pool
        ,E2  = 2*parms0$km                   ##<< total enzyme pool
        ,S1 = 1200                 ##<< N rich substrate
        ,SN1 = 1200/parms0$cnIS1    ##<< N rich substrate N pool
        ,S2 = 170                  ##<< N poor substrate
        ,SN2 = 170/parms0$cnIS2    ##<< N poor substrate N pool
        ,I =  0                    ##<< inorganic pool
)
x <- x0

x0 <- x0Orig <- c( #aE = 0.001*365
        B = 20                     ##<< microbial biomass 
        ,E1  = 1.5*parms0$km                  ##<< total enzyme pool
        ,E2  = 1.5*parms0$km                   ##<< total enzyme pool
        ,S1 = 7000                 ##<< N rich substrate
        ,SN1 = 7000/parms0$cnIS1    ##<< N rich substrate N pool
        ,S2 = 200                  ##<< N poor substrate
        ,SN2 = 200/parms0$cnIS2    ##<< N poor substrate N pool
        ,I =  0                    ##<< inorganic pool
)
x <- x0


derivEezy5 <- function(t,x,parms){
    x <- pmax(unlist(x),1e-16)      # no negative masses
    respMaint <- parms$m * x["B"]
    synE <- parms$aE * x["B"]
    # including maintenance
    E1 <- x["E1"]
    E2 <- x["E2"]
    tvrE1 <- parms$kN1*E1
    tvrE2 <- parms$kN2*E2
    #    
    cnS1 <- x["S1"]/x["SN1"]
    cnS2 <- x["S2"]/x["SN2"]
    limE1 <- E1 / (parms$km1 + E1)
    limE2 <- E2 / (parms$km2 + E2)
    decC1p <- parms$kS1  *  x["S1"]
    decC2p <- parms$kS2  *  x["S2"] 
    decC1 <- decC1p * limE1
    decC2 <- decC2p * limE2
    #decC1 <- parms$kS1  * E1 / (parms$km + E1) * 1000#*  x["S1"]
    #decC2 <- parms$kS2  * E2 / (parms$km + E2) * 1000#*  x["S2"]
    uscE1 <- parms$kN*x["E1"]        # production necessary to balance decay of E1 (uptake steady carbon E1)
    uscE2 <- parms$kN*x["E2"]     
    # average cnE Required for enzyme production according to current concentrations (usc weighted)
    # cnE <- (parms$cnE1*uscE1 + parms$cnE2*uscE2 ) / (uscE1+uscE2)
    tvrERecycling = parms$kNB*(tvrE1+tvrE2)
    uC <- decC1 + decC2 + tvrERecycling
    decN <- decC1/cnS1 + decC2/cnS2 + tvrERecycling/parms$cnE
    plantNUp <- pmin(parms$plantNUp, decN/2)    # plants get at maximum halv of the available N
    uN <-  decN - plantNUp
    # C balance under C limitation
    effE1 <- decC1/uscE1
    effE2 <- decC2/uscE2
    alphaC <- effE1/(effE1 + effE2)
    # average cN reqruired according to enzyme allocation accoording to C efficiencies
    cnE <- parms$cnE #alphaC*parms$cnE1 + (1-alphaC)*parms$cnE2
    #synE <- uC*parms$aE
    uCGrowth <- uC - synE/parms$eps - respMaint    # C for both growth and growth respiration
    if( uCGrowth > 0){
        synB <- parms$eps*uCGrowth
        respSynB <- (1-parms$eps)*uCGrowth
    }else{
        synB <- uCGrowth    # negative
        respSynB <- 0
    }
    synBC <- synB       # used in calculation of alpha
    uNReq <- synE/cnE + synB/parms$cnB
    effE1 <- (decC1/cnS1)/(uscE1/parms$cnE1)
    effE2 <- (decC2/cnS2)/(uscE2/parms$cnE2)
    alphaN <- effE1/(effE1 + effE2)
    if( uN < uNReq ){
        # N limitation
        # average cnE reqruired according to enzyme allocation accoording to N efficiencies
        #cnE <- alpha*parms$cnE1 + (1-alpha)*parms$cnE2
        #synE <- min( uC*parms$aE, uN*cnE )
        uNSynB <- uN - synE/cnE  
        synB <- uNSynB*parms$cnB  # C for both growth and growth respiration
        if( synB > 0){
            respSynB <- (1-parms$eps)/parms$eps*synB
        }else{
            respSynB <- 0
        }
    }
    respSynE <- (1-parms$eps)/parms$eps * synE
    tvrB <- parms$tau*x["B"]
    respTvr <- (1-parms$epsTvr) * tvrB 
    respO <- uC- (synE+synB+respSynE+respSynB+respMaint)
    MmImb <- uN - (synE/parms$cnE + synB/parms$cnB)
    MmTvr <- (1-parms$epsTvr)*tvrB/parms$cnB
    Mm <- MmImb + MmTvr
    #
    pNLim <- (uNReq/uN)^20
    pCLim <- ((uC+respO)/uC)^20
    #alpha <- pCLim*alphaC + pNLim*alphaN / (pCLim + pNLim)
    alpha <- (pCLim*alphaC + pNLim*alphaN) / (pCLim + pNLim)
    if( isTRUE(parms$isAlphaMatch) ){
        synB0 <- max(0, synBC)
        #cnOpt <- (parms$cnE*synE + parms$cnB*synB0)/(synE+synB0)  # optimal biomass ratio
        cnOpt <- (synE+synB0)/(synE/cnE + synB0/parms$cnB)  # optimal biomass ratio
        # calcMatchAlphaEnz
        alpha <- calcMatchAlpha( E=E1+E2, decPotS1=decC1p, decPotS2=decC2p, respMaint=respMaint 
                ,cnOpt=cnOpt
                , cnS1 = cnS1, cnS2=cnS2
                , parms=parms)
    }
    if( isTRUE(parms$isAlphaFix) ){
        alpha <- 0.5
    }        
    #
    # tvr that feeds back to S1 pool, assume that N in SOM for resp (by epsTvr) is mineralized
    tvrC <- (1-parms$kNB)*(tvrE1 +tvrE2) + parms$epsTvr*tvrB
    tvrN <- (1-parms$kNB)*(tvrE1 +tvrE2)/parms$cnE + parms$epsTvr*tvrB/parms$cnB
    #
    dB <- synB - tvrB
    dE1 <- + alpha*synE  - tvrE1
    dE2 <- + (1-alpha)*synE  - tvrE2
    if( isTRUE(parms$isFixedS) ){
        # scenario of fixed substrate
        dS1 <- dS2 <- dSN1 <- dSN2 <- 0
        tvrExC <- 0
        tvrExN <- 0
    }else{ 
        dS2 <- +parms$iS2 -decC2
        dSN2 <- +parms$iS2/parms$cnIS2 -decC2/cnS2
        if( isTRUE(parms$isTvrNil) ){
            # scenario of enzymes and biomass not feeding back to S1
            dS1 <- +parms$iS1 -decC1
            dSN1 <- +parms$iS1/parms$cnIS1 -decC1/cnS1
            tvrExC <- tvrC
            tvrExN <- tvrN
        }else{
            dS1 <- +parms$iS1 +tvrC -decC1
            dSN1 <- +parms$iS1/parms$cnIS1 +tvrN -decC1/cnS1
            tvrExC <- 0
            tvrExN <- 0
        }
    }
    dI <- +Mm
    respB <- respSynE + respSynB + respMaint + respO 
    resp <- respB + respTvr
    #
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dSN1, dS2, dSN2, dI)),names=c("dB","dE1","dE2","dS1","dSN1","dS2","dSN2","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    sqrEps <- sqrt(.Machine$double.eps)
    if( diff( unlist(c(uC=uC, usage=respB+synB+synE )))^2 > sqrEps )  stop("biomass mass balance C error")
    if( diff( unlist(c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + MmImb )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    if( !isTRUE(parms$isFixedS) ){    
        if( diff(unlist(c(dB+dE1+dE2+dS1+dS2+resp+tvrExC,    parms$iS1+parms$iS2 )))^2 > sqrEps )  stop("mass balance C error")
        if( diff(unlist(c((dB)/parms$cnB+(dE1+dE2)/parms$cnE+dSN1+dSN2+dI+tvrExN,    parms$iS1/parms$cnIS1+parms$iS2/parms$cnIS2-plantNUp)))^2 > .Machine$double.eps )  stop("mass balance N error")
    }
    list( resDeriv, c(respO=as.numeric(respO)
        , Mm=as.numeric(Mm), MmImb=as.numeric(MmImb), MmTvr=as.numeric(MmTvr)  
        , alpha=as.numeric(alpha)
        , isLimN=as.numeric(uNReq/uN), uNReq=as.numeric(uNReq), uN=as.numeric(uN)
        , cnS1=as.numeric(cnS1), cnS2=as.numeric(cnS2)
        , limE1=as.numeric(limE1), limE2=as.numeric(limE2)
        , decC1=as.numeric(decC1), decC2=as.numeric(decC2)
        , resp=as.numeric(resp), respB=as.numeric(respB), respTvr=as.numeric(respTvr)
        , tvrB=as.numeric(tvrB)
    ))
}

calcAlphaMatch <- function( cnOpt, B, d1p, d2p, aE, kN1, kN2, kNB, cnS1, cnS2, cnE){
    
    
}

calcS2Steady5 <- function(
    # calculate S2 steady state for given input and size of S1 pool 
    i2
    ,S1
    ,parms
    ,alpha=0.5
){
    with(parms,
       #i2*(S1*aE*alpha*eps*kS1 + aE*alpha*eps*i2 - aE*kN1*km1 + 2*aE*kN1*km2 - eps*kN1*km1*m + 2*eps*kN1*km2*m - kN1*km1*tau + 2*kN1*km2*tau - sqrt(S1^2*aE^2*alpha^2*eps^2*kS1^2 + 2*S1*aE^2*alpha^2*eps^2*i2*kS1 - 2*S1*aE^2*alpha*eps*kN1*kS1*km1 - 2*S1*aE*alpha*eps^2*kN1*kS1*km1*m - 2*S1*aE*alpha*eps*kN1*kS1*km1*tau+ aE^2*alpha^2*eps^2*i2^2 + 2*aE^2*alpha*eps*i2*kN1*km1 + aE^2*kN1^2*km1**2 + 2*aE*alpha*eps^2*i2*kN1*km1*m + 2*aE*alpha*eps*i2*kN1*km1*tau + 2*aE*eps*kN1^2*km1^2*m + 2*aE*kN1^2*km1^2*tau + eps^2*kN1^2*km1^2*m^2 + 2*eps*kN1**2*km1^2*m*tau + kN1^2*km1^2*tau^2))/(kS2*(S1*aE*alpha*eps*kS1 + aE*alpha*eps*i2 - aE*kN1*km1 - eps*kN1*km1*m - kN1*km1*tau - sqrt(S1^2*aE^2*alpha^2*eps**2*kS1^2 + 2*S1*aE^2*alpha^2*eps^2*i2*kS1 - 2*S1*aE^2*alpha*eps*kN1*kS1*km1 - 2*S1*aE*alpha*eps^2*kN1*kS1*km1*m - 2*S1*aE*alpha*eps*kN1*kS1*km1*tau + aE^2*alpha^2*eps^2*i2^2 + 2*aE^2*alpha*eps*i2*kN1*km1 + aE^2*kN1^2*km1^2 + 2*aE*alpha*eps^2*i2*kN1*km1*m + 2*aE*alpha*eps*i2*kN1*km1*tau + 2*aE*eps*kN1^2*km1^2*m + 2*aE*kN1^2*km1^2*tau + eps^2*kN1^2*km1^2*m^2 + 2*eps*kN1^2*km1^2*m*tau + kN1^2*km1^2*tau^2)))            
       i2*(S1*aE*alpha*eps*kS1 + aE*alpha*eps*i2 - aE*kN1*km1 + 2*aE*kN1*km2 - eps*kN1*km1*m + 2*eps*kN1*km2*m - kN1*km1*tau + 2*kN1*km2*tau + sqrt(S1^2*aE^2*alpha^2*eps^2*kS1^2 + 2*S1*aE^2*alpha^2*eps^2*i2*kS1 - 2*S1*aE^2*alpha*eps*kN1*kS1*km1 - 2*S1*aE*alpha*eps^2*kN1*kS1*km1*m - 2*S1*aE*alpha*eps*kN1*kS1*km1*tau+ aE^2*alpha^2*eps^2*i2^2 + 2*aE^2*alpha*eps*i2*kN1*km1 + aE^2*kN1^2*km1^2 + 2*aE*alpha*eps^2*i2*kN1*km1*m + 2*aE*alpha*eps*i2*kN1*km1*tau + 2*aE*eps*kN1^2*km1^2*m + 2*aE*kN1^2*km1^2*tau + eps^2*kN1^2*km1^2*m^2 + 2*eps*kN1^2*km1^2*m*tau + kN1^2*km1^2*tau^2))/(kS2*(S1*aE*alpha*eps*kS1 + aE*alpha*eps*i2 - aE*kN1*km1 - eps*kN1*km1*m - kN1*km1*tau + sqrt(S1^2*aE^2*alpha^2*eps^2*kS1^2 + 2*S1*aE^2*alpha^2*eps^2*i2*kS1 - 2*S1*aE^2*alpha*eps*kN1*kS1*km1 - 2*S1*aE*alpha*eps^2*kN1*kS1*km1*m - 2*S1*aE*alpha*eps*kN1*kS1*km1*tau + aE^2*alpha^2*eps^2*i2^2 + 2*aE^2*alpha*eps*i2*kN1*km1 + aE^2*kN1^2*km1^2 + 2*aE*alpha*eps^2*i2*kN1*km1*m + 2*aE*alpha*eps*i2*kN1*km1*tau + 2*aE*eps*kN1^2*km1^2*m + 2*aE*kN1^2*km1^2*tau + eps^2*kN1^2*km1^2*m^2 + 2*eps*kN1^2*km1^2*m*tau + kN1^2*km1^2*tau^2)))    
       )
}

calcBSteady5 <- function(
        # calculate S2 steady state for given input and size of S1 pool 
        i2
        ,S1
        ,S2 = calcS2Steady5(i2,S1,parms,alpha)
        ,parms
        ,alpha=0.5
){
    with(parms,
            (S1*aE*alpha*eps*kS1 + aE*alpha*eps*i2 - aE*kN1*km1 - eps*kN1*km1*m - kN1*km1*tau + sqrt(S1^2*aE^2*alpha^2*eps^2*kS1^2 + 2*S1*aE^2*alpha^2*eps^2*i2*kS1- 2*S1*aE^2*alpha*eps*kN1*kS1*km1 - 2*S1*aE*alpha*eps^2*kN1*kS1*km1*m - 2*S1*aE*alpha*eps*kN1*kS1*km1*tau + aE^2*alpha^2*eps^2*i2^2 + 2*aE^2*alpha*eps*i2*kN1*km1 + aE^2*kN1^2*km1^2 + 2*aE*alpha*eps^2*i2*kN1*km1*m + 2*aE*alpha*eps*i2*kN1*km1*tau + 2*aE*eps*kN1^2*km1^2*m + 2*aE*kN1^2*km1^2*tau + eps^2*kN1^2*km1^2*m^2 + 2*eps*kN1^2*km1^2*m*tau + kN1^2*km1^2*tau^2))/(2*aE*alpha*(aE + eps*m + tau))    
    )
}

calcESteady5 <- function(
        # calculate S2 steady state for given input and size of S1 pool 
        B
        ,parms
        ,alpha=0.5
){
    with(parms,
            c(alpha, 1-alpha) * aE * B / c(kN1,kN2)
    )
}




.tmp.f <- function(){
    #x0["S2"] <- calcS2Steady5( parms0$iS2, x0["S1"], parms=parms0 )
    #x0["B"] <- calcBSteady5( parms0$iS2, x0["S1"], S2=x0["S2"], parms=parms0 )
    #x0[c("E1","E2")] <- calcESteady5( x0["B"], parms=parms0 )
    
    parmsTvr0 <- within(parms0,{
                isTvrNil <- TRUE
                iS1 <- 160
            })
    parmsInit <- parms0
    #parmsInit <- within(parms0, {isFixedS <- TRUE})
    #parmsInit <- within(parms0, {isAlphaFix <- TRUE})
    #parmsInit <- within(parms0, {isAlphaMatch <- TRUE})
    derivEezy5(0, x0, parmsInit)
    times <- seq(0,2100, length.out=101)
    #times <- seq(0,10000, length.out=101)
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parmsInit))
    #res <- res1f <- as.data.frame(lsoda( x0, times, derivEezy5, parms=within(parms0, useFixedAlloc<-TRUE) ))
    xE <- unlist(tail(res,1))
    plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
    #plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100","E1_1000"))
    #calcS2Steady5( parms0$iS2, xE["S1"], parms=parms0, alpha=xE["alpha"] )
    #calcBSteady5( parms0$iS2, xE["S1"], parms=parms0, alpha=xE["alpha"] )
    #trace(derivEezy5, recover) #untrace(derivEezy5)
    tmp <- derivEezy5(0, xE[1:length(x0)+1], parmsInit)
    
    parmsC2 <- within(parmsInit, { # double as much C in 
                iS2 <- iS2*2
                cnIS2 <- cnIS2*2
                #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
            }) 
    times <- seq(0,10, length.out=101)
    res <- res2S <- as.data.frame(lsoda( xE[1:length(x0)+1], times, derivEezy5, parms=parmsC2))
    xE2 <- unlist(tail(res,1))
    plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
    #trace(derivEezy5, recover)        #untrace(derivEezy5)
    derivEezy5(0, xE2[1:length(x0)+1], parmsC2)
    
    res <- res3S <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivEezy5, parms=parmsInit))
    plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
    xE3 <- tail(res,1)
    
    parmsC0 <- within(parmsInit, { # double as much C in 
                iS2 <- iS2/100
                #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
            }) 
    times <- seq(0,30000, length.out=101)
    res <- res2S <- as.data.frame(lsoda( xE[1:length(x0)+1], times, derivEezy5, parms=parmsC0))
    xE4 <- unlist(tail(res,1))
    plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
    #trace(derivEezy5, recover)        #untrace(derivEezy5)
    derivEezy5(0, xE2[1:length(x0)+1], parmsC2)
    
    
    
 
    data.frame( alphaBest=res$alpha, alphaReal=res$E1/(res$E1+res$E2))
    
    
}




