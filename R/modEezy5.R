#library(deSolve)
# feeding back microbial turnover and enzymes to the S1 pool, 
# now explicitly modelling N because CN-ratio can change
# 

# mg C and N, days
parms0 <- list(
        cnB = 7.16
        ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        ,cnS1 = 5
        ,cnS2 = 100 ##<< N poor substrate
        ,kN = 0.05   ##<< (per day) enzyme turnover
        ,kS1 = 2e-3      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        ,kS2 = 10e-3      ##<< substrate decomposition rate N-poor
        #,aE = 0.05   ##<< C-uptake allocated to enzymes
        ,aE = 0.005   ##<< C-uptake allocated to enzymes
        ,K = 0.3     ##<< enzyme half-saturation constant
        ,m = 0.01    ##<< maintenance respiration rate
        ,tau = 0.012    ##<< biomass turnover rate
        ,eps = 0.5      ##<< carbon use efficiency
        ,iS1 = 0        ##<< here input by enzyme and biomass tvr are modelled explicitely
        ,iS2 = 0.5     
)
parms0 <- within(parms0,{
 K1 <- K2 <- K
 eps1 <- eps2 <- eps
 cnE1 <- cnE2 <- cnE 
        })

parms <- parms0

x0 <- x0Orig <- c(
        B = 10            ##<< microbial biomass 
        ,E1  = 0.01        ##<< total enzyme pool
        ,E2  = 0.01        ##<< total enzyme pool
        ,S1 = 200          ##<< N rich substrate C pool
        ,SN1 = 200/parms0$cnS1          ##<< N rich substrate N pool
        ,S2 = 200         ##<< N poor substrate (cn is that of input and not changing)
        ,I =  0            ##<< inorganic pool
)
x <- x0

derivEezy5 <- function(t,x,parms){
    respMaint <- parms$m * x["B"]
    synE <- parms$aE * x["B"]
    # including maintenance
    E1 <- x["E1"]
    E2 <- x["E2"]
    cnS1 <- x["S1"]/x["SN1"]
    decC1 <- parms$kS1  * E1 / (parms$K + E1) *  x["S1"]
    decC2 <- parms$kS2  * E2 / (parms$K + E2) *  x["S2"]
    #decC1 <- parms$kS1  * E1 / (parms$K + E1) * 1000#*  x["S1"]
    #decC2 <- parms$kS2  * E2 / (parms$K + E2) * 1000#*  x["S2"]
    uscE1 <- parms$kN*x["E1"]        # production necessary to balance decay of E1 (uptake steady carbon E1)
    uscE2 <- parms$kN*x["E2"]     
    # average cnE Required for enzyme production according to current concentrations (usc weighted)
    # cnE <- (parms$cnE1*uscE1 + parms$cnE2*uscE2 ) / (uscE1+uscE2)  
    uC <- decC1 + decC2
    uN <- decC1/cnS1 + decC2/parms$cnS2
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
    uNReq <- synE/cnE + synB/parms$cnB
    effE1 <- (decC1/cnS1)/(uscE1/parms$cnE1)
    effE2 <- (decC2/parms$cnS2)/(uscE2/parms$cnE2)
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
    respO <- uC- (synE+synB+respSynE+respSynB+respMaint)
    Mm <- uN - (synE/parms$cnE + synB/parms$cnB)
    #
    pNLim <- (uNReq/uN)^10
    pCLim <- (uC/(uC+respO))^10
    alpha <- pCLim*alphaC + pNLim*alphaN / (pCLim + pNLim)
    #
    tvrE1 <- parms$kN*x["E1"]
    tvrE2 <- parms$kN*x["E2"]
     #    
    dB <- synB - tvrB
    dE1 <- + alpha*synE  - tvrE1
    dE2 <- + (1-alpha)*synE  - tvrE2
    dS1 <- +parms$iS1 +tvrE1 +tvrE2 +tvrB -decC1
    dSN1 <- +parms$iS1/parms$cnS1 +(tvrE1 +tvrE2)/parms$cnE +tvrB/parms$cnB -decC1/cnS1
    dS2 <- +parms$iS2 -decC2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    #
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dSN1, dS2, dI)),names=c("dB","dE1","dE2","dS1","dSN1","dS2","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    if( diff( unlist(c(uC=uC, usage=resp+synB+synE )))^2 > .Machine$double.eps )  stop("biomass mass balance C error")
    if( diff( unlist(c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + Mm )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    if( diff(unlist(c(dB+dE1+dE2+dS1+dS2+resp,    parms$iS1+parms$iS2 )))^2 > .Machine$double.eps )  stop("mass balance C error")
    if( diff(unlist(c((dB)/parms$cnB+(dE1+dE2)/parms$cnE+dSN1+dS2/parms$cnS2+dI,    parms$iS1/parms$cnS1+parms$iS2/parms$cnS2 )))^2 > .Machine$double.eps )  stop("mass balance N error")
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm), alpha=as.numeric(alpha), isLimN=as.numeric(ifelse(uN < uNReq,1,0)), cnS1=as.numeric(cnS1) ))
}


.tmp.f <- function(){
    times <- seq(0,5000, length.out=101)
    derivEezy5(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    parmsMsmall <- within(parms0, m <- 1e-3)   # small maintenance
    res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parmsM0))
    
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parms0))
    res <- res1b <- as.data.frame(lsoda( x0, times, derivEezy5, parms=within(parms0,useAlpha0<-TRUE)))
    
    res <- res2 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parmsMsmall))
    
    #res <- lsoda( x0, times, derivEezy, parms=parms0)
    
    #trace(derivEezy4, recover)  # untrace(derivEezy4)
    derivEezy5(0, tail(res,1)[2:7], parmsM0)
    derivEezy5(0, tail(res,1)[2:7], parms0)
    
    parms1 <- within(parms0, kS2<-0 )
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parms1))

    
    
    res$B1000 <- res$B/1000
    res$B100 <- res$B/100
    res$E1_10 <- res$E1/10
    res$E2_10 <- res$E2/10
    res$S1f <- res$S1/x0["S1"]
    res$S2f <- res$S2/x0["S2"]
    res$cnS1_10 <- res$cnS1/10
    cls <- c("B100","E1_10","E2_10","respO","Mm","S1f","S2f","cnS1_10")
    #cls <- c("B1000","E1_10","E2_10","Mm")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:10, col=1:10)
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
 
    data.frame( alphaBest=res$alpha, alphaReal=res$E1/(res$E1+res$E2))
}




