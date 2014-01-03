#library(deSolve)
# extension to have really two enzyme pools, only synthesis is split 
# and better partitioning of flux
# no biosynthesis if more C required for enzyme and maintenance than taken up ->

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
        #,DOM = 0           ##<< dissolved organic matter
        ,E1  = 0.01  ##<< total enzyme pool
        ,E2  = 0.01  ##<< total enzyme pool
        ,S1 = 500       ##<< N rich substrate
        ,S2 = 500         ##<< N poor substrate
        ,I =  0            ##<< inorganic pool
)
x <- x0



parms0 <- list(
        cnB = 7.16
        ,cnE = 7.16     #tom 3
        ,cnS1 = 5
        ,kN = 0.05   ##<< (per day) enzyme turnover
        ,kS = 1     ##<< substrate decomposition rate
        ,aE = 0.05   ##<< C-uptake allocated to enzymes
        ,K = 0.3    ##<< enzyme half-saturation constant
        ,m = 0.01    ##<< maintenance respiration rate
        ,tau = 0.012    ##<< biomass turnover rate
        ,eps = 0.5      ##<< carbon use efficiency
)
parms0 <- within(parms0,{
 kS1 <- kS2 <- kS
 K1 <- K2 <- K
 eps1 <- eps2 <- eps
        })

parms <- parms0
x <- x0

derivEezy2 <- function(t,x,parms){
    E <- x["E1"]+x["E2"]
    alpha0 <- x["E1"]/E  
    respMaint <- parms$m * x["B"]
    # including maintenance 
    alpha  <- alphaM <- 1/2 - sqrt(
                    -E^2*(E + 2*parms$K)*(2*parms$eps*parms$kS - respMaint)*(-2*E*parms$eps*parms$kS + E*respMaint + 2*parms$K*respMaint)
                ) /
                (2*E^2*(2*parms$eps*parms$kS - respMaint))
    if( !is.finite(alpha) | (alpha <0) | (alpha > 1) ) 
        alpha  <- 1/2 * (2 * parms$eps * parms$cnS1 * E-parms$cnB * parms$K-parms$cnB * E+(4 * (parms$eps * parms$cnS1 * E)^2 -8 * parms$eps * parms$cnS1 * E * parms$cnB * parms$K-4 * parms$eps * parms$cnS1 * E^2 * parms$cnB+(parms$cnB * parms$K)^2 +2 * parms$cnB^2 * parms$K * E+(parms$cnB * E)^2+8 * (parms$eps * parms$cnS1)^2 * parms$K * E)^(1/2) ) /  E/(2 * parms$eps * parms$cnS1-parms$cnB) 
    #alpha2 <- 1/2 + sqrt(-E^2*(E + 2*parms$K)*(2*parms$eps*parms$kS - respMaint)*(-2*E*parms$eps*parms$kS + E*respMaint + 2*parms$K*respMaint))/(2*E^2*(2*parms$eps*parms$kS - respMaint))    
    E1 <- x["E1"]
    E2 <- x["E2"]
    decC1 <- parms$kS *  E1 / (parms$K + E1)
    decC2 <- parms$kS *  E2 / (parms$K + E2)
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    minUCE <- min(uC, uN*parms$cnE)  # C for enzyme biosynthesis, if necessary take microbial material for associated resp
    synE <- parms$aE * minUCE
    respSynE <- (1-parms$eps)/parms$eps * synE
    minUCM <- min( uC-synE-respSynE-respMaint, (uN-synE/parms$cnE)*parms$cnB/parms$eps ) # C flux available for biomass synthesis
    if( minUCM > 0){    
        synB <- parms$eps * minUCM  
        respSynB <- (1-parms$eps)/parms$eps * synB
        edecay <- 0
    }else{
        synB <- minUCM      # negative: endogeneous decay, releases N
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
    resDeriv <- structure(c( dB, dE1, dE2, dS1, dS2, dI),names=c("dB","dE1","dE2","dS1","dS2","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm)), alpha=as.numeric(alpha) )
}


.tmp.f <- function(){
    times <- seq(0,300, length.out=101)
    derivEezy2(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parmsM0))
    #res <- lsoda( x0, times, derivEezy, parms=parms0)
    
    #trace(derivEezy2, recover)  # untrace(derivEezy2)
    derivEezy2(0, tail(res,1)[2:7], parms0)
    
    parms1 <- within(parms0, kS2<-0 )
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parms1))

    
    
    res$B1000 <- res$B/1000
    res$E10 <- res$E1/10
    res$E210 <- res$E2/10
    cls <- c("B1000","E10","E210","respO","Mm")
#cls <- c("B","E","respO","Mm","S1","S2")
    matplot( res[,1], res[,cls], type="l")
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
 
    data.frame( alphaBest=res$alpha, alphaReal=res$E1/(res$E1+res$E2))
}




