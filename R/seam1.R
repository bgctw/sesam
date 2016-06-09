#library(deSolve)
# Refactoring of Easy5 with terminology corresponding to SoilPaper14 

# gC/m2 and gN/m2, /yr

derivSeam1 <- function(t,x,parms){
    x <- pmax(unlist(x),1e-16)      # no negative masses
    #
    ER <- x["ER"]
    EL <- x["EL"]
    cnR <- x["R"]/x["RN"]
    cnL <- x["L"]/x["LN"]
    cnE <- parms$cnE #alphaC*parms$cnER + (1-alphaC)*parms$cnEL
    cnB <- parms$cnB
    #
    synE <- parms$aE * x["B"]       # total enzyme production
    rM <- parms$m * x["B"]          # maintenance respiration
    tvrB <- parms$tau*x["B"]        # microbial turnover
    tvrER <- parms$kNR*ER
    tvrEL <- parms$kNL*EL
    #    
    decRp <- parms$kR * x["R"]
    decLp <- parms$kL * x["L"] 
    limER <- ER / (parms$kmR + ER)
    limEL <- EL / (parms$kmL + EL)
    decR <- decRp * limER
    decL <- decLp * limEL
    #
    tvrERecycling = parms$kNB*(tvrER+tvrEL)
    uC <- decR + decL + tvrERecycling
    CsynBC <- uC - synE/parms$eps - rM
    #
    # Nitrogen balance
    decN <- decR/cnR + decL/cnL + tvrERecycling/parms$cnE
    plantNUp <- pmin(parms$plantNUp, decN/2)    # plants get at maximum halv of the available N
    uN <-  decN - plantNUp
    NsynBN <- uN - synE/cnE
    CsynBN <- (cnB * NsynBN)/parms$eps      # C available for biomass growth and growth respiration
    #
    isLimN <- CsynBN <  CsynBC 
    CsynB = if( isLimN ) CsynBN else CsynBC
    #
    # average cN reqruired according to enzyme allocation accoording to C efficiencies
    #synE <- uC*parms$aE
    #uCGrowth <- uC - synE/parms$eps - rM    # C for both growth and growth respiration
    if( CsynB > 0){
        synB <- parms$eps*CsynB
        rG <- (1-parms$eps)*CsynB
    }else{
        synB <- CsynB    # negative
        rG <- 0
    }
    # imbalance fluxes
    respSynE <- (1-parms$eps)/parms$eps * synE
    respO <- uC- (synE+respSynE+synB+rG+rM)
    MmImb <- uN - (synE/parms$cnE + synB/parms$cnB)
    respTvr <- (1-parms$epsTvr) * tvrB 
    MmTvr <- respTvr/parms$cnB
    Mm <- MmImb + MmTvr
    #
    # Revenue strategy
    revLC <- decL / (parms$kNL) / (parms$kmL + EL)
    revRC <- decR / (parms$kNR) / (parms$kmR + ER)
    revLN <- revLC * cnE/cnL
    revRN <- revRC * cnE/cnR
    alphaC <- revRC / (revLC + revRC)  
    alphaN <- revRN / (revLN + revRN)  
    #
    delta <- 200
    wCLim = min( .Machine$double.xmax, (CsynBN/CsynBC)^delta )      # min to avoid +Inf
    wNLim = min( .Machine$double.xmax, (parms$eps*CsynBC/cnB / NsynBN)^delta )
    alpha <- (wCLim*alphaC + wNLim*alphaN) / (wCLim + wNLim) 
    
    if( isTRUE(parms$isAlphaMatch) ){
        synB0 <- max(0, synB)
        #cnOpt <- (parms$cnE*synE + parms$cnB*synB0)/(synE+synB0)  # optimal biomass ratio
        cnOpt <- (synE+synB0)/(synE/cnE + synB0/parms$cnB)  # optimal biomass ratio
        # calcMatchAlphaEnz
        alpha <- calcMatchAlphaSeam( E=ER+EL, decPotR=decRp, decPotL=decLp, respMaint=rM 
                ,cnOpt=cnOpt
                , cnR = cnR, cnL=cnL
                , parms=parms)
    }
    if( isTRUE(parms$isAlphaFix) ){
        alpha <- 0.5
    }        
    #
    # tvr that feeds back to R pool, assume that N in SOM for resp (by epsTvr) is mineralized
    tvrC <- +parms$epsTvr*tvrB  +(1-parms$kNB)*(tvrER +tvrEL) 
    tvrN <- +parms$epsTvr*tvrB/parms$cnB  +(1-parms$kNB)*(tvrER +tvrEL)/parms$cnE 
    tvrExC <- 0     # fluxes leaving the system
    tvrExN <- 0
    #
    dB <- synB - tvrB
    dER <- + alpha*synE  - tvrER
    dEL <- + (1-alpha)*synE  - tvrEL
    dL <- -decL +parms$iL 
    dLN <- -decL/cnL  +parms$iL/parms$cnIL 
    dR <- -decR +parms$iR +tvrC 
    dRN <- -decR/cnR +parms$iR/parms$cnIR +tvrN 
    #
    if( isTRUE(parms$isFixedS) ){
        # scenario of fixed substrate
        dR <- dL <- dRN <- dLN <- 0
    }else{ 
        if( isTRUE(parms$isTvrNil) ){
            # scenario of enzymes and biomass not feeding back to R
            dR <- +parms$iR -decR
            dRN <- +parms$iR/parms$cnIR -decR/cnR
            tvrExC <- tvrC
            tvrExN <- tvrN
        }
    }
    dI <- +Mm
    respB <- respSynE + rG + rM + respO 
    resp <- respB + respTvr
    #
    resDeriv <-structure(as.numeric(c( dB, dER, dEL, dR, dRN, dL, dLN, dI)),names=c("dB","dER","dEL","dR","dRN","dL","dLN","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    sqrEps <- sqrt(.Machine$double.eps)
    # parms$iL - (decL + dL)
    # parms$iR + +tvrC -(decR + dR)
    # 
    if( diff( unlist(c(uC=uC, usage=respB+synB+synE )))^2 > sqrEps )  stop("biomass mass balance C error")
    if( diff( unlist(c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + MmImb )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    if( !isTRUE(parms$isFixedS) ){    
        if( diff(unlist(c(dB+dER+dEL+dR+dL+resp+tvrExC,    parms$iR+parms$iL )))^2 > sqrEps )  stop("mass balance C error")
        if( diff(unlist(c((dB)/parms$cnB+(dER+dEL)/parms$cnE+dRN+dLN+dI+tvrExN,    parms$iR/parms$cnIR+parms$iL/parms$cnIL-plantNUp)))^2 > .Machine$double.eps )  stop("mass balance N error")
    }
    if( isTRUE(parms$isRecover) ) recover()    
    list( resDeriv, c(respO=as.numeric(respO)
        , Mm=as.numeric(Mm), MmImb=as.numeric(MmImb), MmTvr=as.numeric(MmTvr)  
        , alpha=as.numeric(alpha)
        , alphaC=as.numeric(alphaC), alphaN=as.numeric(alphaN)
        , cnR=as.numeric(cnR), cnL=as.numeric(cnL)
        , limER=as.numeric(limER), limEL=as.numeric(limEL)
        , decR=as.numeric(decR), decL=as.numeric(decL)
        , resp=as.numeric(resp), respB=as.numeric(respB), respTvr=as.numeric(respTvr)
        , tvrB=as.numeric(tvrB)
        , revRC=as.numeric(revRC), revLC=as.numeric(revLC), revRN=as.numeric(revRN), revLN=as.numeric(revLN)
        , pCsyn = as.numeric(CsynBC / CsynBN), CsynReq=as.numeric(CsynBN), Csyn=as.numeric(CsynBC)
        , pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) ), NsynReq=as.numeric(CsynBC/cnB), Nsyn=as.numeric(NsynBN)
    #    wCLim = (CsynBN/CsynBC)^delta
    #    wNLim = (parms$eps*CsynBC/cnB / NsynBN)^delta
))
}








plotResSeam1 <- function(res, legendPos="topleft"
        , cls = c("B","ER","EL","respO","Mm")
        , xlab="time (yr)"
        , ylab="gC/m2 , gN/m2 , % ,/yr"
){
    #res$B100 <- res$B/100
    res$ER_10 <- res$ER*10
    res$EL_10 <- res$EL*10
    res$ER_1000 <- res$ER*1000
    res$EL_1000 <- res$EL*1000
    res$alpha100 <- res$alpha*100
    res$B100 <- res$B * 100
    res$B10 <- res$B * 10
    res$Rr <- res$R / res$R[1]  * 100 #max(res$R, na.rm=TRUE)
    res$Lr <- res$L / res$L[1]  * 100 #max(res$L, na.rm=TRUE)
    #cls <- c("B100","ER_10","EL_10","respO","Mm")
    #cls <- c("B100","ER_10","EL_10","respO","Mm","eff")
    #cls <- c("B1000","ER_10","EL_10","Mm")
    #cls <- c("B","E","respO","Mm","R","L")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:20, col=1:20, xlab=xlab, ylab="")
    mtext( ylab, 2, line=2.5, las=0)
    legend(legendPos, inset=c(0.01,0.01), legend=cls, lty=1:20, col=1:20)
}

