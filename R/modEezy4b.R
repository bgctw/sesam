#library(deSolve)
# overall biomass ratio 
# compromize between C limited and N limited enzyme allocation quota 
# enzyme production proportional to biomass instead of uptake
# no changing S1 and S2

# mg C and N, days
# moving to gC/m2 yr


x0 <- x0Orig <- c(
        B = 30             ##<< microbial biomass 
        ,E1  = 3        ##<< total enzyme pool
        ,E2  = 3        ##<< total enzyme pool
        ,S1 = 2000    ##<< N rich substrate
        ,S2 = 260    ##<< N poor substrate
)
x <- x0



parms0 <- list(
        cnB = 7.16
        #,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        ,cnE = 7.16
        ,cnS1 = 7
        ,cnS2 = 70 ##<< N poor substrate
        #,kN = 0.05   ##<< (per day) enzyme turnover
        ,kN= 0.05*365  ##<< /yr enzyme turnover 5% turning over each day
        #,kS1 = 0.2      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        #,kS2 = 1      ##<< substrate decomposition rate N-poor
        #,kS1 = 5e-3      ##<< substrate decomposition rate N-rich (here simulate large N stock)
        #,kS2 = 10e-3     ##<< substrate decomposition rate N-poor
        #,kS2 = 5e-2     ##<< substrate decomposition rate N-poor
        #,aE = 0.05   ##<< C-uptake allocated to enzymes
        #,kS1 = 1/(5*365)        ##<< 5 years 
        #,kS2 = 1/(0.5*365)        ##<< 1/2 years 
        ,kS1 = 1/(5)        ##<< 5 years 
        ,kS2 = 1/(0.5)        ##<< 1/2 years 
        ,aE = 0.01*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
        #,km = 0.3     ##<< enzyme half-saturation constant
        ,km = 2     ##<< enzyme half-saturation constant
        ,m = 0.01*365    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = 0.012*365    ##<< biomass turnover rate
        ,eps = 0.5      ##<< carbon use efficiency
        ,iS1 = 220           #xE["B"]*(parms$aE+parms$tau)
        ,iS2 = 300         # g/m2 input per year (half NPP)
        ,useFixedAlloc=FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
)
parms0 <- within(parms0,{
 km1 <- km2 <- km
 eps1 <- eps2 <- eps
 kN1 <- kN2 <- kN
 cnE1 <- cnE2 <- cnE 
        })

parms <- parms0
x <- x0

derivEezy4b <- function(t,x,parms){
    x <- pmax(unlist(x),1e-16)      # no negative masses
    respMaint <- parms$m * x["B"]
    synE <-  parms$aE * x["B"]
    # including maintenance
    E1 <- x["E1"]
    E2 <- x["E2"]
    #decC1 <- parms$kS1S  * E1 / (parms$km + E1) #* 1000#*  x["S1"]
    #decC2 <- parms$kS2S  * E2 / (parms$km + E2) #* 1000#*  x["S2"]
    limE1 <- E1 / (parms$km1 + E1)
    limE2 <- E2 / (parms$km2 + E2)
    decC1 <- parms$kS1  * limE1 *  x["S1"]
    decC2 <- parms$kS2  * limE2 *  x["S2"]
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
    dS1 <- +parms$iS1 -decC1
    dS2 <- +parms$iS2 -decC2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    effE1=as.numeric(ifelse(uN < uNReq,effE1C,effE1N))
    effE2=as.numeric(ifelse(uN < uNReq,effE2C,effE2N))
    #
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dS2)),names=c("dB","dE1","dE2","dS1","dS2"))
    #resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dS2, dI)),names=c("dB","dE1","dE2","dS1","dS2","dI"))
    #resDeriv <-structure(as.numeric(c( dB, dE1, dE2)),names=c("dB","dE1","dE2"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    if( diff( unlist(c(uC=uC, usage=resp+synB+synE )))^2 > .Machine$double.eps )  stop("biomass mass balance C error")
    if( diff( unlist(c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + Mm )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    if( diff(unlist(c(dB+dE1+dE2+dS1+dS2+tvrE1+tvrE2+resp+tvrB,    parms$iS1+parms$iS2 )))^2 > .Machine$double.eps )  stop("mass balance C error")
    if( diff(unlist(c((dB+tvrB)/parms$cnB+(dE1+dE2+tvrE1+tvrE2)/parms$cnE+dS1/parms$cnS1+dS2/parms$cnS2+dI,    parms$iS1/parms$cnS1+parms$iS2/parms$cnS2 )))^2 > .Machine$double.eps )  stop("mass balance N error")
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm), alpha=as.numeric(alpha)
        , isLimN=as.numeric(ifelse(uN < uNReq,1,0)), respSynE=as.numeric(respSynE), respSynB=as.numeric(respSynB), respMaint=as.numeric(respMaint) 
        , effE1=effE1  , effE2=effE2, eff=as.numeric(alpha*effE1+(1-alpha)*effE2)
        , limE1=as.numeric(limE1), limE2=as.numeric(limE2)
))
}


.tmp.f <- function(){
    times <- seq(0,500, length.out=101)
    derivEezy4b(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    parmsMsmall <- within(parms0, m <- 1e-3)   # small maintenance
    #res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy4b, parms=parmsM0))
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy4b, parms=parms0))
    
    res$B1000 <- res$B/1000
    res$B100 <- res$B/100
    res$B10 <- res$B/10
    res$E1_10 <- res$E1/10
    res$E2_10 <- res$E2/10
    res$S1_10 <- res$S1/10
    res$S2_10 <- res$S2/10
    cls <- c("B10","E1","E2","respO","Mm","S1_10","S2_10")
    #cls <- c("B1000","E1_10","E2_10","Mm")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:10, col=1:10, ylim=c(0,70), xlab="time (d)", ylab="")
    legend("topright", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
 
    aE <- aEs[3]
    data.frame( alphaBest=res$alpha, alphaReal=res$E1/(res$E1+res$E2))
    with( as.list(tail(res,1)),   c(alpha*aE*B ,  E1*parms$kN )
    )
    
    tmp <- with( c(parmsAe, as.list(tail(res,1))),   list(
                    B=B
                    ,a=aE, tau1=kN, tau2=kN
                    ,d1=kS1
                    ,d2=kS2
                    ,m1=km, m2=km
                    ,cn1=cnS1, cn2=cnS2, cnE=cnE
    ))
    alphaC <- with(tmp,
            (2*B*a*d1 + d1*m2*tau2 + d2*m1*tau1 - sqrt(4*B^2*a^2*d1*d2 + 4*B*a*d1*d2*m1*tau1 + 4*B*a*d1*d2*m2*tau2 + d1^2*m2^2*tau2^2 + 2*d1*d2*m1*m2*tau1*tau2 + d2^2*m1^2*tau1^2))/(2*B*a*d1 - 2*B*a*d2)
    )
    alphaN <- with(tmp,
            (2*B*a*cn2*d1 + cn1*d2*m1*tau1 + cn2*d1*m2*tau2 - sqrt(4*B^2*a^2*cn1*cn2*d1*d2 + 4*B*a*cn1*cn2*d1*d2*m1*tau1 + 4*B*a*cn1*cn2*d1*d2*m2*tau2 + cn1^2*d2^2*m1^2*tau1^2 + 2*cn1*cn2*d1*d2*m1*m2*tau1*tau2 + cn2^2*d1^2*m2^2*tau2^2))/(-2*B*a*cn1*d2 + 2*B*a*cn2*d1)
    )
    
    
}

.tmp.f.Steady <- function(){
    x0sL <- estSteadyState4bc(parms0, d1p=parms0$kS1*x0["S1"]*0.8, d2p=parms0$kS2*x0["S2"]*0.8)
    x0s <- x0
    x0s[names(x0sL$x0)] <-  x0sL$x0   
    #trace(derivEezy4bc, recover)   #untrace(derivEezy4bc)
    derivEezy4b(0, x0s, parms0)
    times <- seq(0,10, length.out=101)
    res <- res1S <- as.data.frame(lsoda( x0s, times, derivEezy4b, parms=parms0))
    xE <- unlist(tail(res,1))
    plotRes(res, "topright", cls = c("B","respO","Mm","S1r","S2r","alpha100"))
    #abline(h=0.5, col="lightgray", lty="dashed")
    
    parmsC2 <- within(parms0, { # double as much C in 
               iS2 <- iS2*2
               cnS2 <- cnS2*2
            }) 
    res <- res2S <- as.data.frame(lsoda( xE[1:length(x0)+1], times, derivEezy4b, parms=parmsC2))
    plotRes(res, "topright", cls = c("B","respO","Mm","S1r","S2r","alpha100"))
    xE2 <- unlist(tail(res,1))
    
    res <- res3S <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivEezy4b, parms=parms0))
    plotRes(res, "topright", cls = c("B","respO","Mm","S1r","S2r","alpha100"))
    xE3 <- tail(res,1)
    
}



plotRes <- function(res, legendPos="topleft"
    , cls = c("B","E1","E2","respO","Mm")
    , xlab="time (yr)"
    , ylab="gC/m2 , gN/m2 , % ,/yr"
){
    #res$B100 <- res$B/100
    res$E1_10 <- res$E1*10
    res$E2_10 <- res$E2*10
    res$E1_1000 <- res$E1*1000
    res$E2_1000 <- res$E2*1000
    res$alpha100 <- res$alpha*100
    res$B100 <- res$B * 100
    res$B10 <- res$B * 10
    res$S1r <- res$S1 / res$S1[1]  * 100 #max(res$S1, na.rm=TRUE)
    res$S2r <- res$S2 / res$S2[1]  * 100 #max(res$S2, na.rm=TRUE)
    #cls <- c("B100","E1_10","E2_10","respO","Mm")
    #cls <- c("B100","E1_10","E2_10","respO","Mm","eff")
    #cls <- c("B1000","E1_10","E2_10","Mm")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:20, col=1:20, xlab=xlab, ylab="")
    mtext( ylab, 2, line=2.5, las=0)
    legend(legendPos, inset=c(0.01,0.01), legend=cls, lty=1:20, col=1:20)
}


.tmp.faE <- function(){
    # test effect of aE on achieved biomass
    aEs <- seq(0.001, 0.02, length.out=20)
    aEi <- aEs[1]
    resL <- lapply( aEs, function(aEi){
                parmsAe <- within(parms, aE <- aEi)
                res <- as.data.frame(lsoda( x0, times, derivEezy4b, parms=parmsAe))
            })
    parmsAe <- within(parms, aE <- 0.003)
    res <- resL[[3]]
    #res <- resL[[4]]
    #res <- resL[[2]]       # less overflow but less biomass
    tail(res,1)
    plotRes( res )
    
    B <- sapply(resL, function(res){ tail(res,1)[["B"]]})
    eff <- sapply(resL, function(res){ tail(res,1)[["eff"]]})
    alpha <- sapply(resL, function(res){ tail(res,1)[["alpha"]]})
    plot( B ~ aEs)
}






