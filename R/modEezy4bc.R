#library(deSolve)
# same as 4b but prescribing fixed potential influx (no changing S1 and S2)

.tmp.f <- function(){
x0 <- x0Orig <- c(
        B = 10             ##<< microbial biomass 
        ,E1  = 0.01        ##<< total enzyme pool
        ,E2  = 0.01        ##<< total enzyme pool
)
x <- x0

parms0 <- list(
        cnB = 7.16
        #,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        ,cnE = 7.16
        #,cnS1 = 5
        ,cnS1 = 8
        ,cnS2 = 100 ##<< N poor substrate
        ,kN = 0.05   ##<< (per day) enzyme turnover
        #,kS1 = 0.2      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        #,kS2 = 1      ##<< substrate decomposition rate N-poor
        ,kS1 = 5e-3      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        ,kS2 = 5e-2     ##<< substrate decomposition rate N-poor
        #,aE = 0.05   ##<< C-uptake allocated to enzymes
        ,aE = 0.01   ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,km = 0.3     ##<< enzyme half-saturation constant
        ,m = 0.01    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = 0.012    ##<< biomass turnover rate
        ,eps = 0.5      ##<< carbon use efficiency
        #,iS1 = 1
        #,iS2 = 1
        ,useFixedAlloc=FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
)
parms0 <- within(parms0,{
            kS1S <- kS1*500
            kS2S <- kS2*100
            km1 <- km2 <- km
            kN1 <- kN2 <- kN
         eps1 <- eps2 <- eps
         cnE1 <- cnE2 <- cnE 
           })

parms <- parms0
x <- x0
}

derivEezy4bc <- function(t,x,parms){
    x <- pmax(unlist(x),1e-16)      # no negative masses
    respMaint <- parms$m * x["B"]
    synE <-  parms$aE * x["B"]
    # including maintenance
    E1 <- x["E1"]
    E2 <- x["E2"]
    decC1 <- parms$kS1S  * E1 / (parms$km + E1) #* 1000#*  x["S1"]
    decC2 <- parms$kS2S  * E2 / (parms$km + E2) #* 1000#*  x["S2"]
    uscE1 <- parms$kN*x["E1"]        # production necessary to balance decay of E1 (uptake steady carbon E1)
    uscE2 <- parms$kN*x["E2"]     
    # average cnE Required for enzyme production according to current concentrations (usc weighted)
    # cnE <- (parms$cnE1*uscE1 + parms$cnE2*uscE2 ) / (uscE1+uscE2)  
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1 + decC2/parms$cnS2
    # C balance under C limitation
    effE1C <- decC1/uscE1
    effE2C <- decC2/uscE2
    alphaC <- effE1C/(effE1C + effE2C)
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
    uNReq <- synE/cnE + synB/parms$cnB
    effE1N <- (decC1/parms$cnS1)/(uscE1/parms$cnE1)
    effE2N <- (decC2/parms$cnS2)/(uscE2/parms$cnE2)
    alphaN <- effE1N/(effE1N + effE2N)
    if( uN < uNReq ){
        # N limitation
        # average cnE reqruired according to enzyme allocation accoording to N efficiencies
        #cnE <- alpha*parms$cnE1 + (1-alpha)*parms$cnE2
        #synE <- min( uC*parms$aE, uN*cnE, parms$aE*x["B"] )
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
    respO <- uC- (synE+synB+respSynE+respSynB+respMaint)
    Mm <- uN - (synE/parms$cnE + synB/parms$cnB)
    #
    pNLim <- (uNReq/uN)^100
    pCLim <- ((uC+respO)/uC)^100
    alpha <- if( isTRUE(parms$useFixedAlloc)) 0.5 else (pCLim*alphaC + pNLim*alphaN) / (pCLim + pNLim)
     #    
    tvrE1 <- parms$kN*x["E1"]
    tvrE2 <- parms$kN*x["E2"]
    #
    dB <- synB - tvrB
    dE1 <- + alpha*synE  - tvrE1
    dE2 <- + (1-alpha)*synE  - tvrE2
    dS1 <- -decC1  #+parms$iS1 
    dS2 <- -decC2  #+parms$iS2 
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    effE1=as.numeric(ifelse(uN < uNReq,effE1C,effE1N))
    effE2=as.numeric(ifelse(uN < uNReq,effE2C,effE2N))
    #
    #resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dS2, dI)),names=c("dB","dE1","dE2","dS1","dS2","dI"))
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2)),names=c("dB","dE1","dE2"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    if( diff( unlist(c(uC=uC, usage=resp+synB+synE )))^2 > .Machine$double.eps )  stop("biomass mass balance C error")
    if( diff( unlist(c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + Mm )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    #if( diff(unlist(c(dB+dE1+dE2+dS1+dS2+tvrE1+tvrE2+resp+tvrB,    parms$iS1+parms$iS2 )))^2 > .Machine$double.eps )  stop("mass balance C error")
    if( diff(unlist(c(dB+dE1+dE2+dS1+dS2+tvrE1+tvrE2+resp+tvrB,    0 )))^2 > .Machine$double.eps )  stop("mass balance C error")
    #if( diff(unlist(c((dB+tvrB)/parms$cnB+(dE1+dE2+tvrE1+tvrE2)/parms$cnE+dS1/parms$cnS1+dS2/parms$cnS2+dI,    parms$iS1/parms$cnS1+parms$iS2/parms$cnS2 )))^2 > .Machine$double.eps )  stop("mass balance N error")
    if( diff(unlist(c((dB+tvrB)/parms$cnB+(dE1+dE2+tvrE1+tvrE2)/parms$cnE+dS1/parms$cnS1+dS2/parms$cnS2+dI,    0 )))^2 > .Machine$double.eps )  stop("mass balance N error")
#recover()    
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm), alpha=as.numeric(alpha), alphaC=as.numeric(alphaC), alphaN=as.numeric(alphaN)
        , isLimN=as.numeric(ifelse(uN < uNReq,1,0)), respSynE=as.numeric(respSynE), respSynB=as.numeric(respSynB), respMaint=as.numeric(respMaint) 
        , effE1=effE1  , effE2=effE2, eff=as.numeric(alpha*effE1+(1-alpha)*effE2)
))
}

alphaSteadyC4b <- function(B, d1p, d2p, parms=parms0){
    with(parms,
            (2*B*aE*d1p + d1p*kN2*km2 + d2p*kN1*km1 - sqrt(4*B^2*aE^2*d1p*d2p + 4*B*aE*d1p*d2p*kN1*km1 + 4*B*aE*d1p*d2p*kN2*km2 + d1p^2*kN2^2*km2^2 + 2*d1p*d2p*kN1*kN2*km1*km2 + d2p^2*kN1^2*km1^2))/(2*B*aE*d1p - 2*B*aE*d2p)
    )
}

alphaSteadyN4b <- function(B, d1p, d2p,parms=parms0){
    with(parms,
            -(2*B*aE*cnS2*d1p + cnS1*d2p*kN1*km1 + cnS2*d1p*kN2*km2 - sqrt(4*B^2*aE^2*cnS1*cnS2*d1p*d2p + 4*B*aE*cnS1*cnS2*d1p*d2p*kN1*km1 + 4*B*aE*cnS1*cnS2*d1p*d2p*kN2*km2 + cnS1^2*d2p^2*kN1^2*km1^2 + 2*cnS1*cnS2*d1p*d2p*kN1*kN2*km1*km2 + cnS2^2*d1p^2*kN2^2*km2^2))/(2*B*aE*cnS1*d2p - 2*B*aE*cnS2*d1p)
    )
}

.fdBC <- function(B,aEOpt=parms$aE, d1p,d2p, parms){
    with( within(parms, aE <- aEOpt),
    B*aE*d1p*(2*B*aE*d1p + d1p*kN2*km2 + d2p*kN1*km1 - sqrt(4*B^2*aE^2*d1p*d2p + 4*B*aE*d1p*d2p*kN1*km1 + 4*B*aE*d1p*d2p*kN2*km2 + d1p^2*kN2^2*km2^2 + 2*d1p*d2p*kN1*kN2*km1*km2 + d2p^2*kN1^2*km1^2))/(kN1*(2*B*aE*d1p - 2*B*aE*d2p)*(B*aE*(2*B*aE*d1p + d1p*kN2*km2 + d2p*kN1*km1 - sqrt(4*B^2*aE^2*d1p*d2p + 4*B*aE*d1p*d2p*kN1*km1 + 4*B*aE*d1p*d2p*kN2*km2 + d1p^2*kN2^2*km2^2 + 2*d1p*d2p*kN1*kN2*km1*km2 + d2p^2*kN1^2*km1^2))/(kN1*(2*B*aE*d1p - 2*B*aE*d2p)) + km1)) + B*aE*d2p*(1 - (2*B*aE*d1p + d1p*kN2*km2 + d2p*kN1*km1 - sqrt(4*B^2*aE^2*d1p*d2p +4*B*aE*d1p*d2p*kN1*km1 + 4*B*aE*d1p*d2p*kN2*km2 + d1p^2*kN2^2*km2^2 + 2*d1p*d2p*kN1*kN2*km1*km2 + d2p^2*kN1^2*km1^2))/(2*B*aE*d1p - 2*B*aE*d2p))/(kN2*(B*aE*(1 - (2*B*aE*d1p + d1p*kN2*km2 + d2p*kN1*km1 - sqrt(4*B^2*aE^2*d1p*d2p + 4*B*aE*d1p*d2p*kN1*km1 + 4*B*aE*d1p*d2p*kN2*km2 + d1p^2*kN2^2*km2^2 + 2*d1p*d2p*kN1*kN2*km1*km2 + d2p^2*kN1^2*km1^2))/(2*B*aE*d1p - 2*B*aE*d2p))/kN2 + km2))- B*aE/eps - B*m - B*tau/eps 
    )
}
.fdBCOptB <- function( B,...){
    .fdBC( B, ... )^2    
}
.fdBCOptAE <- function(  aE  ,... ){
    res <- optimize( .fdBCOptB, c(1,1000), aEOpt=aE, ...)
    res$minimum
}
attr(.fdBC, "ex") <- function(){
    Bs <- seq(1,200,by=1)
    plot( .fdBC(Bs, aEOpt=0.01, d1p=d1p, d2p=d2p, parms=parms) ~ Bs, type="l" )
    abline(h=0, col="grey", lty="dashed")
    abline(v=136)
    #res <- optimize( .fdBCOptB, c(1,1000), aEOpt=0.01, d1p=d1p, d2p=d2p, parms=parms, tol=1e-13)
    res <- optimize( .fdBCOptB, c(1,1000), aEOpt=0.01, d1p=d1p, d2p=d2p, parms=parms)
    #.fdBCOptAE( 0.01, d1p=d1p, d2p=d2p, parms=parms )
    res <- optimize( .fdBCOptAE, c(0.001,1), d1p=d1p, d2p=d2p, parms=parms, maximum = TRUE)
    abline(v=res$objective, col="blue")
    lines( .fdBC(Bs, aEOpt=res$maximum, d1p=d1p, d2p=d2p, parms=parms) ~ Bs, col="blue" )
    c(aEOpt=res$maximum)
}



.fdBN <- function(B,aEOpt=parms$aE, d1p,d2p,parms){
    with( parms,
    B*aE*d2p*(1 - (-2*B*aE*cnS2*d1p - cnS1*d2p*kN1*km1 - cnS2*d1p*kN2*km2 + sqrt(4*B^2*aE^2*cnS1*cnS2*d1p*d2p + 4*B*aE*cnS1*cnS2*d1p*d2p*kN1*km1 + 4*B*aE*cnS1*cnS2*d1p*d2p*kN2*km2 + cnS1^2*d2p^2*kN1^2*km1^2 + 2*cnS1*cnS2*d1p*d2p*kN1*kN2*km1*km2 + cnS2^2*d1p^2*kN2^2*km2^2))/(2*B*aE*cnS1*d2p - 2*B*aE*cnS2*d1p))/(cnS2*kN2*(B*aE*(1 - (-2*B*aE*cnS2*d1p - cnS1*d2p*kN1*km1 - cnS2*d1p*kN2*km2 + sqrt(4*B^2*aE^2*cnS1*cnS2*d1p*d2p + 4*B*aE*cnS1*cnS2*d1p*d2p*kN1*km1 + 4*B*aE*cnS1*cnS2*d1p*d2p*kN2*km2 + cnS1^2*d2p^2*kN1^2*km1^2 + 2*cnS1*cnS2*d1p*d2p*kN1*kN2*km1*km2 + cnS2^2*d1p^2*kN2^2*km2^2))/(2*B*aE*cnS1*d2p - 2*B*aE*cnS2*d1p))/kN2 + km2)) + B*aE*d1p*(-2*B*aE*cnS2*d1p - cnS1*d2p*kN1*km1 - cnS2*d1p*kN2*km2+ sqrt(4*B^2*aE^2*cnS1*cnS2*d1p*d2p + 4*B*aE*cnS1*cnS2*d1p*d2p*kN1*km1 + 4*B*aE*cnS1*cnS2*d1p*d2p*kN2*km2 + cnS1^2*d2p^2*kN1^2*km1^2 + 2*cnS1*cnS2*d1p*d2p*kN1*kN2*km1*km2 + cnS2^2*d1p^2*kN2^2*km2^2))/(cnS1*kN1*(2*B*aE*cnS1*d2p - 2*B*aE*cnS2*d1p)*(B*aE*(-2*B*aE*cnS2*d1p - cnS1*d2p*kN1*km1 - cnS2*d1p*kN2*km2 +sqrt(4*B^2*aE^2*cnS1*cnS2*d1p*d2p + 4*B*aE*cnS1*cnS2*d1p*d2p*kN1*km1 + 4*B*aE*cnS1*cnS2*d1p*d2p*kN2*km2 + cnS1^2*d2p^2*kN1^2*km1^2 + 2*cnS1*cnS2*d1p*d2p*kN1*kN2*km1*km2 + cnS2^2*d1p^2*kN2^2*km2^2))/(kN1*(2*B*aE*cnS1*d2p - 2*B*aE*cnS2*d1p)) + km1)) - B*aE/cnE - B*tau/cnB
    )
}
.fdBNOptB <- function( B,...){
    .fdBN( B, ... )^2    
}

estSteadyState4bc <- function(
        ### estimate steady biomass for given decomposition rates and parameters
        parms                 
        ,d1p=parms$kS1S     ##<< potential decomposition flux from S1 (nutrient rich)
        ,d2p=parms$kS2S     ##<< potential decomposition flux from S2 (nutrient poor)
        ,aE=parms$aE        ##<< enzyme allocation constant
){
    res <- optimize( .fdBCOptB, c(1,1000), aEOpt=aE, d1p=d1p, d2p=d2p, parms=parms)
    BC <- res$minimum
    res <- optimize( .fdBNOptB, c(1,1000), aEOpt=aE, d1p=d1p, d2p=d2p, parms=parms)
    BN <- res$minimum
    isLimN <- (BN <= BC)
    B <- if(isLimN) BN else BC
    alpha <- if(isLimN) alphaSteadyN4b(BN, d1p=d1p, d2p=d2p, parms) else alphaSteadyC4b(BC,d1p=d1p, d2p=d2p, parms)
    allocE <- parms$aE * B / parms$kN1
    E1 <- alpha * allocE
    E2 <- (1-alpha) * allocE
    list(
            x0 = c(B=B, E1=as.numeric(E1), E2=as.numeric(E2))
            ,isLimN = isLimN
            ,alpha = alpha
            )
}






.tmp.f <- function(){
    times <- seq(0,500, length.out=101)
    derivEezy4bc(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    parmsMsmall <- within(parms0, m <- 1e-3)   # small maintenance
    #res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy4b, parms=parmsM0))
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy4bc, parms=parms0))
    
    res$B1000 <- res$B/1000
    res$B100 <- res$B/100
    res$E1_10 <- res$E1/10
    res$E2_10 <- res$E2/10
    cls <- c("B100","E1_10","E2_10","respO","Mm")
    #cls <- c("B1000","E1_10","E2_10","Mm")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:20, col=1:20)
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:20, col=1:20)
    
    xE <- tail(res,1)
    x0sL <- estSteadyState4bc(parms0)
    x0s <-  x0sL$x0   
    #trace(derivEezy4bc, recover)   #untrace(derivEezy4bc)
    derivEezy4bc(0, xE, parms0)
    derivEezy4bc(0, x0s, parms0)
    res <- res1S <- as.data.frame(lsoda( x0s, times, derivEezy4bc, parms=parms0))
    
}




