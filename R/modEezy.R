#library(deSolve)

# original Eezy model
# both overflow and mineralization?

# mg C and N, days
tmp.f.depr <- function(){
    initVars <- list(
            B = c(1,0.14)            ##<< microbial biomass 
            ,DOM = c(0,0)           ##<< dissolved organic matter
            ,E  = c( 0.02, 2.8e-3)  ##<< total enzyme pool
            ,S1 = c(500, 100)       ##<< N rich substrate
            ,S2 = c(500, 0)         ##<< N poor substrate
            ,I = c(0, 0)            ##<< inorganic pool
    )
    xInit <- do.call( rbind, initVars )  #, names=c("C","N") )
    colnames(xInit) <- c("C","N")
}

.tmp.f <- function(){
	x0 <- x0Orig <- c(
	        B = 10          ##<< microbial biomass 
	        #,DOM = 0       ##<< dissolved organic matter
	        ,E  = 0.02      ##<< total enzyme pool
	        ,S1 = 500       ##<< N rich substrate
	        ,S2 = 500       ##<< N poor substrate
	        ,I =  0         ##<< inorganic pool
	)
	x <- x0
	
	
	
	parms0 <- list(
	        cnB = 7.16
	        ,cnE = 7.16     #tom 3
	        ,cnS1 = 5
	        ,kN = 0.05   ##<< (per day) enzyme turnover
	        ,kS = 1      ##<< substrate decomposition rate
	        ,aE = 0.05   ##<< C-uptake allocated to enzymes
	        ,K = 0.3     ##<< enzyme half-saturation constant
	        ,m = 0.01    ##<< maintenance respiration rate
	        ,tau = 0.012    ##<< biomass turnover rate
	        ,eps = 0.5      ##<< carbon use efficiency
	)
	parms0 <- within(parms0,{
	 kS1 <- kS2 <- kS*500
	 K1 <- K2 <- K
	 eps1 <- eps2 <- eps
	        })
	
	parms <- parms0
	x <- x0
}

derivEezyOrig <- function(t,x,parms){
    # only one enzme pool E_Total is simulated. at each time this split into E1 and E2 by alpha
    alpha  <- 1/2 * (2 * parms$eps * parms$cnS1 * x["E"]-parms$cnB * parms$K-parms$cnB * x["E"]+(4 * (parms$eps * parms$cnS1 * x["E"])^2 -8 * parms$eps * parms$cnS1 * x["E"] * parms$cnB * parms$K-4 * parms$eps * parms$cnS1 * x["E"]^2 * parms$cnB+(parms$cnB * parms$K)^2 +2 * parms$cnB^2 * parms$K * x["E"]+(parms$cnB * x["E"])^2+8 * (parms$eps * parms$cnS1)^2 * parms$K * x["E"])^(1/2) ) /  x["E"]/(2 * parms$eps * parms$cnS1-parms$cnB)  
    E1 <- alpha * x["E"]
    E2 <- (1-alpha) * x["E"]
    decC1 <- parms$kS *  E1 / (parms$K + E1)
    decC2 <- parms$kS *  E2 / (parms$K + E2)
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    minUC <- min(uC, uN*parms$cnE) # C for enzyme biosynthesis, if necessary take microbial material for associated resp
    minUCM <- min(uC, uN*parms$cnB) # C for biomass biosynthesis
    synE <- parms$aE * minUC
    respSynE <- (1-parms$eps)/parms$eps * synE
    respMaint <- parms$m * x["B"]
    synB <- parms$eps * (1-parms$aE/parms$eps) * (minUCM - respMaint)
    #synB2 <- min( parms$eps*(uC - synE - respSynE -respMaint), (uN-synE/parms$cnE)*parms$cnB   )
    respSynB <- (1-parms$eps)/parms$eps * synB
    tvrB <- parms$tau*x["B"]
    respO <- uC - (synE+synB+respSynE+respSynB+respMaint)
    Mm <- uN - (synE/parms$cnE + synB/parms$cnB)
            
    dB <- +synB -tvrB
    dE <- +synE  - parms$kN*x["E"]
    dS1 <- -decC1
    dS2 <- -decC2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    
    # c(decC=decC, usage=resp+synB+synE )  # check the same
    # c(decN=decN, usage=synE/pr$cnE + synB/pr$cnB + Mm )  # check the same
    
    list( c( dB, dE, dS1, dS2, dI), c(respO=as.numeric(respO), Mm=as.numeric(Mm)) )
}

derivEezy <- function(t,x,parms){
    # changed formulation of synB to take minimum of C or N limitation
    alpha = 0.5
    alpha  <- 1/2 * (2 * parms$eps * parms$cnS1 * x["E"]-parms$cnB * parms$K-parms$cnB * x["E"]+(4 * (parms$eps * parms$cnS1 * x["E"])^2 -8 * parms$eps * parms$cnS1 * x["E"] * parms$cnB * parms$K-4 * parms$eps * parms$cnS1 * x["E"]^2 * parms$cnB+(parms$cnB * parms$K)^2 +2 * parms$cnB^2 * parms$K * x["E"]+(parms$cnB * x["E"])^2+8 * (parms$eps * parms$cnS1)^2 * parms$K * x["E"])^(1/2) ) /  x["E"]/(2 * parms$eps * parms$cnS1-parms$cnB)
    if( alpha < 0 ) alpha = 0
    if( alpha > 1 ) alpha = 1
    E1 <- alpha * x["E"]
    E2 <- (1-alpha) * x["E"]
    decC1 <- parms$kS1 *  E1 / (parms$K + E1)
    decC2 <- parms$kS2 *  E2 / (parms$K + E2)
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    minUC <- min(uC, uN*parms$cnE) # C for enzyme biosynthesis, if necessary take microbial material for associated resp
    synE <- parms$aE * minUC              
    respSynE <- (1-parms$eps)/parms$eps * synE
    respMaint <- parms$m * x["B"]
    synB <- min( parms$eps*(uC - synE - respSynE -respMaint), (uN-synE/parms$cnE)*parms$cnB   )
    respSynB <- (1-parms$eps)/parms$eps * synB
    tvrB <- parms$tau*x["B"]
    respO <- uC - (synE+synB+respSynE+respSynB+respMaint)
    Mm <- uN - (synE/parms$cnE + synB/parms$cnB)
    
    dB <- +synB -tvrB
    dE <- +synE  - parms$kN*x["E"]
    dS1 <- -decC1
    dS2 <- -decC2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    
    # c(decC=decC, usage=resp+synB+synE )  # check the same
    # c(decN=decN, usage=synE/pr$cnE + synB/pr$cnB + Mm )  # check the same
    
    list( c( dB, dE, dS1, dS2, dI), c(respO=as.numeric(respO), Mm=as.numeric(Mm), decC1=as.numeric(decC1), decC2=as.numeric(decC2) ) )
}



.tmp.f <- function(){
    times <- seq(0,1000, length.out=101)
    derivEezyOrig(0, x0, parms0)
    res <- res0 <- as.data.frame(lsoda( x0, times, derivEezyOrig, parms=parms0))
    #res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy, parms=parms0))
    #res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy_2, parms=parms0))
    #res <- lsoda( x0, times, derivEezy, parms=parms0)
    
    #trace(derivEezyOrig, recover)   # untrace(derivEezyOrig)
    derivEezyOrig(0, tail(res,1)[2:6], parms0)

    
    # with changed synB formulation
    times <- seq(0,1000, length.out=101)
    derivEezy(0, x0, parms0)
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy, parms=parms0))
    
    
    .tmp.C1only <- function(){
        parms1 <- within(parms, kS2<-0 )
        res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy_2, parms=parms1))
    }
        
    res$B1000 <- res$B/1000
    res$E10 <- res$E/10
    cls <- c("B1000","E10","respO","Mm")
#cls <- c("B","E","respO","Mm","S1","S2")
    matplot( res[,1], res[,cls], type="l")
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
    
}

.tmp.fCNGraph <- function(){
    cnS1s <- seq(5,25,by=.2)
    
    .tmp.f <- function(){
        # no decomposition of S2        
        parmsNoS2 <- within(parms0, kS2 <- 0)
        cnAll <- cnS1s
        resl <- lapply(cnS1s, function(cnS1i){
                    parms <- within(parmsNoS2, cnS1 <- cnS1i)
                    res <- as.data.frame(lsoda( x0, times, derivEezy, parms=parms))
                    tail(res,1)
                })
    }
    
    
    #cnS1i <- cnS1s[1]
    #cnS1i <- cnS1s[length(cnS1s)]
    cnAll <- 1000/(500/cnS1s)
    resl <- lapply(cnS1s, function(cnS1i){
                parms <- within(parms0, cnS1 <- cnS1i)
                res <- as.data.frame(lsoda( x0, times, derivEezy, parms=parms))
                tail(res,1)
            })
    .tmp.f <- function(){
        #trace(derivEezy, recover)   # untrace(derivEezy)
        derivEezy( 0, tail(res,1)[2:6], parms)
    }
    
    resE1 <- do.call( rbind, resl)
    resE <- cbind( cnS1 = cnAll, resE1 )
    resE$B1000 <- resE$B/1000
    resE$E10 <- resE$E/10
    cls <- c("B1000","respO","Mm")
    #cls <- c("B1000","E10","respO","Mm")
    matplot( resE[,1], resE[,cls], type="l")
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
    
    
    matplot( resE)
}




