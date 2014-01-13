#library(deSolve)
# extension to have really two enzyme pools, only synthesis is split 
# partitioning according to efficiencies and limitations split for enzyme and for biomass
# mg C and N, days
tmp.f.depr <- function(){
    initVars <- list(
            B = c(1,0.14)            ##<< microbial biomass 
            ,DOM = c(0,0)           ##<< dissolved organic matter
            ,E1  = c( 0.01, 1.4e-3)  ##<< total enzyme pool
            ,E2  = c( 0.01, 1.4e-3)  ##<< total enzyme pool
            ,S1 = c(500, 100)       ##<< N rich substrate
            ,S2 = c(500, 0)         ##<< N poor substrate
            ,I = c(0, 0)            ##<< inorganic pool
    )
    xInit <- do.call( rbind, initVars )  #, names=c("C","N") )
    colnames(xInit) <- c("C","N")
}


x0 <- x0Orig <- c(
        B = 10            ##<< microbial biomass 
        ,E1  = 0.01        ##<< total enzyme pool
        ,E2  = 0.01        ##<< total enzyme pool
        ,S1 = 2000          ##<< N rich substrate
        ,S2 = 1000         ##<< N poor substrate
        ,I =  0            ##<< inorganic pool
)
x <- x0



parms0 <- list(
        cnB = 7.16
        ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        ,cnS1 = 5
        ,cnS2 = 100 ##<< N poor substrate
        ,kN = 0.05   ##<< (per day) enzyme turnover
        ,kS1 = 5e-3      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        ,kS2 = 10e-3      ##<< substrate decomposition rate N-poor
        ,aE = 0.05   ##<< C-uptake allocated to enzymes
        ,K = 0.3     ##<< enzyme half-saturation constant
        ,m = 0.01    ##<< maintenance respiration rate
        ,tau = 0.012    ##<< biomass turnover rate
        ,eps = 0.5      ##<< carbon use efficiency
        ,useAlpha0 = FALSE  ##<< flag to indicate equation without maintenance
)
parms0 <- within(parms0,{
 K1 <- K2 <- K
 eps1 <- eps2 <- eps
 cnE1 <- cnE2 <- cnE 
        })

parms <- parms0
x <- x0

derivEezy3 <- function(t,x,parms){
    respMaint <- parms$m * x["B"]
    # including maintenance
    E1 <- x["E1"]
    E2 <- x["E2"]
    decC1 <- parms$kS1  * E1 / (parms$K + E1) * 1000#*  x["S1"]
    decC2 <- parms$kS2  * E2 / (parms$K + E2) * 1000#*  x["S2"]
    uscE1 <- parms$kN*x["E1"]        # production necessary to balance decay of E1 (uptake steady carbon E1)
    uscE2 <- parms$kN*x["E2"]        
    # average cnE Required for current enzymes (usc weighted)
    cnE <- (parms$cnE1*uscE1 + parms$cnE2*uscE2 ) / (uscE1+uscE2)  
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1 + decC2/parms$cnS2
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    
    parms$aE*
    if( uN*cnE >= uC){
        # C limited
        # invest into enzymes proportional to their C return effiQciency
        minUCE <- uC
        effE1 <- decC1/uscE1
        effE2 <- decC2/uscE2
    }else{
        # N limited
        # invest into enzymes proportional to their N return efficiency
        minUCE <- uN*cnE
        effE1 <- (decC1/parms$cnS1)/(uscE1/parms$cnE1)
        effE2 <- (decC2/parms$cnS2)/(uscE2/parms$cnE2)
    }
    alpha <- effE1/(effE1 + effE2)
    synE <- parms$aE * minUCE
    respSynE <- (1-parms$eps)/parms$eps * synE
    ucmC <- uC-synE-respSynE-respMaint
    ucmN <- (uN-synE/parms$cnE)*parms$cnB/parms$eps
    isLimN <- if(ucmC < ucmN) 0 else 1
    minUCM <- min( ucmC, ucmN ) # C flux available for biomass synthesis
    if( minUCM > 0){    
        synB <- parms$eps * minUCM  
        respSynB <- (1-parms$eps)/parms$eps * synB
    }else{
        synB <- minUCM     # N from decayed biomass becomes available
        respSynB <- 0
    }
    tvrB <- parms$tau*x["B"]
    respO <- uC - (synE+synB+respSynE+respSynB+respMaint)
    Mm <- uN - (synE/parms$cnE + synB/parms$cnB)
         
    dB <- synB - tvrB
    dE1 <- + alpha*synE  - parms$kN*x["E1"]
    dE2 <- + (1-alpha)*synE  - parms$kN*x["E2"]
    dS1 <- -decC1
    dS2 <- -decC2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    
    # c(uC=uC, usage=resp+synB+synE )  # check the same
    # c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + Mm )  # check the same
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dS2, dI)),names=c("dB","dE1","dE2","dS1","dS2","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm)), alpha=as.numeric(alpha), isLimEN=ifelse(uN*cnE < uC,1,0), isLimN=isLimN )
}


.tmp.f <- function(){
    times <- seq(0,1000, length.out=101)
    derivEezy3(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    parmsMsmall <- within(parms0, m <- 1e-3)   # small maintenance
    res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy3, parms=parmsM0))
    
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy3, parms=parms0))
    res <- res1b <- as.data.frame(lsoda( x0, times, derivEezy3, parms=within(parms0,useAlpha0<-TRUE)))
    
    res <- res2 <- as.data.frame(lsoda( x0, times, derivEezy3, parms=parmsMsmall))
    
    #res <- lsoda( x0, times, derivEezy, parms=parms0)
    
    #trace(derivEezy3, recover)  # untrace(derivEezy3)
    derivEezy3(0, tail(res,1)[2:7], parmsM0)
    derivEezy3(0, tail(res,1)[2:7], parms0)
    
    parms1 <- within(parms0, kS2<-0 )
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy3, parms=parms1))

    
    
    res$B1000 <- res$B/1000
    res$E1_10 <- res$E1/10
    res$E2_10 <- res$E2/10
    cls <- c("B1000","E1_10","E2_10","respO","Mm")
    #cls <- c("B1000","E1_10","E2_10","Mm")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l")
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
 
    data.frame( alphaBest=res$alpha, alphaReal=res$E1/(res$E1+res$E2))
}




