#library(deSolve)
# extension to have really two enzyme pools, only synthesis is split 
# and better partitioning of flux
# no biosynthesis if more C required for enzyme and maintenance than taken up ->

# problems with switching alpha between two values

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


.tmp.f <- function(){
x0 <- x0Orig <- c(
        B = 10       ##<< microbial biomass 
        #,DOM = 0    ##<< dissolved organic matter
        ,E1  = 0.01  ##<< total enzyme pool
        ,E2  = 0.01  ##<< total enzyme pool
        ,S1 = 500    ##<< N rich substrate
        ,S2 = 500    ##<< N poor substrate
        ,I =  0      ##<< inorganic pool
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
        ,useAlpha0 = FALSE  ##<< flag to indicate equation without maintenance
)
parms0 <- within(parms0,{
 kS1 <- kS2 <- kS
 K1 <- K2 <- K
 eps1 <- eps2 <- eps
        })

parms <- parms0
x <- x0
}

derivEezy2 <- function(t,x,parms){
    E <- x["E1"]+x["E2"]
    alpha0 <- x["E1"]/E  
    respMaint <- parms$m * x["B"]
    # including maintenance
    .tmp.f <- function(){
        .sqrt <-  suppressWarnings(sqrt(-E^2*(E + 2*parms$K)*(2*parms$kS - respMaint)*(-2*E*parms$kS + E*respMaint + 2*parms$K*respMaint))/(2*E^2*(2*parms$kS - respMaint)))
        alpha <- if( !is.finite(.sqrt) ){   
            alpha0                
        }else{
            if( .sqrt > 1/2 ) .sqrt <- 1/2
            if( alpha0 > 1/2 ) 1/2+.sqrt else 1/2 - .sqrt
        }
        alpha <- (E*parms$cnB*parms$cnS1 - parms$K*parms$cnB*parms$cnS1 + sqrt(E^2*parms$cnB^2*parms$cnS1^2 + 2*E*parms$K*parms$cnB^2*parms$cnS1^2 + parms$K^2*parms$cnB^2*parms$cnS1^2))/(2*E*parms$cnB*parms$cnS1)
    }
    #alpha <- (E*parms$cnB - 2*E*parms$cnS1*parms$eps + parms$K*parms$cnB - sqrt(E^2*parms$cnB^2 - 4*E^2*parms$cnB*parms$cnS1*parms$eps + 4*E^2*parms$cnS1^2*parms$eps^2 + 2*E*parms$K*parms$cnB^2 - 8*E*parms$K*parms$cnB*parms$cnS1*parms$eps + 8*E*parms$K*parms$cnS1^2*parms$eps^2 + parms$K^2*parms$cnB^2))/(2*E*parms$cnB - 4*E*parms$cnS1*parms$eps)
    alpha <- alphaEq <- -(-E^2*parms$cnB*parms$kS + 2*E^2*parms$cnS1*parms$eps*parms$kS - E^2*parms$cnS1*parms$eps*respMaint - E*parms$K*parms$cnB*parms$kS + sqrt(E^2*parms$cnB^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$kS - 4*E^2*parms$cnS1*parms$eps*parms$kS + 2*E^2*parms$cnS1*parms$eps*respMaint)
    #alpha2 <- -(-E^2*parms$cnB*parms$kS + 2*E^2*parms$cnS1*parms$eps*parms$kS - E^2*parms$cnS1*parms$eps*respMaint - E*parms$K*parms$cnB*parms$kS - sqrt(E^2*parms$cnB^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$kS - 4*E^2*parms$cnS1*parms$eps*parms$kS + 2*E^2*parms$cnS1*parms$eps*respMaint)
    alpha0 <- 1/2 * (2 * parms$eps * parms$cnS1 * E-parms$cnB * parms$K-parms$cnB * E+(4 * (parms$eps * parms$cnS1 * E)^2 -8 * parms$eps * parms$cnS1 * E * parms$cnB * parms$K-4 * parms$eps * parms$cnS1 * E^2 * parms$cnB+(parms$cnB * parms$K)^2 +2 * parms$cnB^2 * parms$K * E+(parms$cnB * E)^2+8 * (parms$eps * parms$cnS1)^2 * parms$K * E)^(1/2) ) /  E/(2 * parms$eps * parms$cnS1-parms$cnB) # take the no-maintenance variant
    if( alpha < 0 ) alpha = 0
    if( alpha > 1 ) alpha = 1
    #alpha <- -(-E^2*parms$cnB*parms$kS - 2*E^2*parms$cnS1*parms$eps*parms$kS - E^2*parms$cnS1*parms$eps*respMaint - E*parms$K*parms$cnB*parms$kS + sqrt(E^2*parms$cnB^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$kS - 4*E^2*parms$cnS1*parms$eps*parms$kS + 2*E^2*parms$cnS1*parms$eps*respMaint)    
    #x[c("B","E1","E2","S1","S2")] <- pmax(0,x[c("B","E1","E2","S1","S2")])
    #alpha2 <- 1/2 + sqrt(-E^2*(E + 2*parms$K)*(2*parms$eps*parms$kS - respMaint)*(-2*E*parms$eps*parms$kS + E*respMaint + 2*parms$K*respMaint))/(2*E^2*(2*parms$eps*parms$kS - respMaint))    
    E1 <- x["E1"]
    E2 <- x["E2"]
    decC1 <- parms$kS *  E1 / (parms$K + E1)
    decC2 <- parms$kS *  E2 / (parms$K + E2)
    parms$eps * (decC1+decC2)/(decC1/parms$cnS1)
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    minUCE <- min(uC, uN*parms$cnE)  # C for enzyme biosynthesis, if necessary take microbial material for associated resp
    isLimEN <- if(uN*parms$cnE < uC) 0 else 1
    synE <- parms$aE * minUCE
    respSynE <- (1-parms$eps)/parms$eps * synE
    ucmC <- uC-synE-respSynE-respMaint
    ucmN <- (uN-synE/parms$cnE)*parms$cnB/parms$eps
    isLimN <- if(ucmC < ucmN) 0 else 1
    minUCM <- min( ucmC, ucmN ) # C flux available for biomass synthesis
    if( minUCM > 0){    
        synB <- parms$eps * minUCM  
        respSynB <- (1-parms$eps)/parms$eps * synB
        edecay <- 0
    }else{
        #print("negative C balance; recover"); recover()
        synB <- minUCM      # negative: endogeneous decay, releases N
        respSynB <- 0
    }
    tvrB <- parms$tau*x["B"]
    respO <- uC - (synE+synB+respSynE+respSynB+respMaint)
    Mm <- uN - (synE/parms$cnE + synB/parms$cnB)
         
    if( synB < tvrB | parms$useAlpha0) alpha <- alpha0  # provide CN/ratio for enzyme production before maintenance
    dB <- synB - tvrB
    dE1 <- + alpha*synE  - parms$kN*x["E1"]
    dE2 <- + (1-alpha)*synE  - parms$kN*x["E2"]
    dS1 <- -decC1
    dS2 <- -decC2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    
    if( abs(diff( cBal <- c(uC=uC, usage=resp+synB+synE ))) > 1e-6 ) stop("C-balance error")  # check the same
    if( abs(diff( cBal <- c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + Mm ) )) > 1e-6) stop("N-balance error") # check the same
    
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dS2, dI)),names=c("dB","dE1","dE2","dS1","dS2","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    #, alpha2 = as.numeric(alpha2)
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm)), alpha=as.numeric(alpha), isLimEN=isLimEN, isLimN=isLimN )
}


.tmp.f <- function(){
    times <- seq(0,250, length.out=101)
    derivEezy2(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    parmsMsmall <- within(parms0, m <- 1e-3)   # small maintenance
    res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parmsM0))
    
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parms0,atol = 1e-4))
    res <- res1b <- as.data.frame(lsoda( x0, times, derivEezy2, parms=within(parms0,useAlpha0<-TRUE)))
    
    res <- res2 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parmsMsmall))
    
    #res <- lsoda( x0, times, derivEezy, parms=parms0)
    
    #trace(derivEezy2, recover)  # untrace(derivEezy2)
    derivEezy2(0, tail(res,1)[2:7], parms0)
    derivEezy2(0, tail(res,1)[2:7], parmsM0)

    
    parms1 <- within(parms0, kS2<-0 )
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parms1))

    parmsMBig <- within(parms0, m <- 0.04 )
    res <- res3 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parmsMBig,atol = 1e-4))
    res <- res30 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=within(parmsMBig,useAlpha0<-TRUE)) )
    
    
    
    res$B1000 <- res$B/1000
    res$E1_10 <- res$E1/10
    res$E2_10 <- res$E2/10
    cls <- c("B1000","E1_10","E2_10","respO","Mm")
    #cls <- c("B1000","E1_10","E2_10","respO","Mm","alpha")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:10, col=1:10)
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
 
    data.frame( alphaBest=res$alpha, alphaReal=res$E1/(res$E1+res$E2))
}




