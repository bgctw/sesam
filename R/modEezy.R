#library(deSolve)

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


x0 <- x0Orig <- c(
        B = 10            ##<< microbial biomass 
        #,DOM = 0           ##<< dissolved organic matter
        ,E  = 0.02  ##<< total enzyme pool
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
derivEezy <- function(t,x,parms){
    alpha = 0.5
    alpha  <- 1/2 * (2 * parms$eps * parms$cnS1 * x["E"]-parms$cnB * parms$K-parms$cnB * x["E"]+(4 * (parms$eps * parms$cnS1 * x["E"])^2 -8 * parms$eps * parms$cnS1 * x["E"] * parms$cnB * parms$K-4 * parms$eps * parms$cnS1 * x["E"]^2 * parms$cnB+(parms$cnB * parms$K)^2 +2 * parms$cnB^2 * parms$K * x["E"]+(parms$cnB * x["E"])^2+8 * (parms$eps * parms$cnS1)^2 * parms$K * x["E"])^(1/2) ) /  x["E"]/(2 * parms$eps * parms$cnS1-parms$cnB)  
    E1 <- alpha * x["E"]
    E2 <- (1-alpha) * x["E"]
    decC1 <- parms$kS *  E1 / (parms$K + E1)
    decC2 <- parms$kS *  E2 / (parms$K + E2)
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    minUC <- min(uC, uN*parms$cnE)       # C for enzyme biosynthesis, if necessary take microbial material for associated resp
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
    
    list( c( dB, dE, dS1, dS2, dI), c(respO=as.numeric(respO), Mm=as.numeric(Mm)) )
}


derivEezy2 <- function(t,x,parms){
    alpha = 0.5
    alpha  <- 1/2 * (2 * parms$eps * parms$cnS1 * x["E"]-parms$cnB * parms$K-parms$cnB * x["E"]+(4 * (parms$eps * parms$cnS1 * x["E"])^2 -8 * parms$eps * parms$cnS1 * x["E"] * parms$cnB * parms$K-4 * parms$eps * parms$cnS1 * x["E"]^2 * parms$cnB+(parms$cnB * parms$K)^2 +2 * parms$cnB^2 * parms$K * x["E"]+(parms$cnB * x["E"])^2+8 * (parms$eps * parms$cnS1)^2 * parms$K * x["E"])^(1/2) ) /  x["E"]/(2 * parms$eps * parms$cnS1-parms$cnB)  
    E1 <- alpha * x["E"]
    E2 <- (1-alpha) * x["E"]
    decC1 <- parms$kS1 *  E1 / (parms$K1 + E1)
    decC2 <- parms$kS2 *  E2 / (parms$K2 + E2)
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    minUC <- min(uC, uN*parms$cnE)       # C for enzyme biosynthesis, if necessary take microbial material for associated resp
    synE <- parms$aE * minUC
    pUC <- minUC/uC
    pU1 <- decC1 / uC 
    pU2 <- decC2 / uC
    epsA <- (pU1*parms$eps1 + pU2*parms$eps2)   # weighted average of both eps
    epsRA <- (pU1*(1-parms$eps1)/parms$eps1 + pU2*(1-parms$eps2)/parms$eps2) # weighted average of fraction assimilated going to resp
    respSynE <-  epsRA * synE
    respMaint <- parms$m * x["B"]
    synB <- min( epsA*(uC - synE - respSynE -respMaint), (uN-synE/parms$cnE)*parms$cnB   )
    respSynB <- epsRA * synB
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


times <- seq(0,1000, length.out=101)
res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parms0))
#res <- lsoda( x0, times, derivEezy, parms=parms0)


parms1 <- within(parms, kS2<-0 )
res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy2, parms=parms1))



res$B1000 <- res$B/1000
res$E10 <- res$E/10
cls <- c("B1000","E10","respO","Mm")
#cls <- c("B","E","respO","Mm","S1","S2")
matplot( res[,1], res[,cls], type="l")
legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)

