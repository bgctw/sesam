#library(deSolve)
# model corresponding to Seam2 with enzyme levels computed by quasi steady state

# gC/m2 and gN/m2, /yr

derivSesam2 <- function(t,x,parms){
    x <- pmax(unlist(x),1e-16)      # no negative masses
    #
	# compute steady state enzyme levels
	decRp <- parms$kR * x["R"]
	decLp <- parms$kL * x["L"] 
	cnR <- x["R"]/x["RN"]
	cnL <- x["L"]/x["LN"]
	cnE <- parms$cnE #alphaC*parms$cnER + (1-alphaC)*parms$cnEL
	cnB <- parms$cnB
	alphaC <- computeAllocationPartitioning( dR=decRp, dL=decLp, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	alphaN <- computeAllocationPartitioning( dR=decRp/cnR, dL=decLp/cnL, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	ETot <- parms$aE * x["B"] / parms$kN
	rM <- parms$m * x["B"]          # maintenance respiration
	tvrB <- parms$tau*x["B"]        # microbial turnover
	synE <- parms$aE * x["B"]       # total enzyme production per microbial biomass
	#
	# declare variables that will be computed/overidden in computeAlphaFluxes
	# else <<- will override bindings at global space
	ER <- EL <- tvrER <- tvrEL <- limER <- limEL <- decR <- decL <- tvrERecycling <- uC <- decC <-  
			CsynBC <- decN <- plantNUp <- PhiU <- immoPot <- uNSubstrate <- uNPot <-
			NsynBN <- CsynBN <- NsynBNSubstrate <- CsynBNSubstrate <- isLimN <- 
			isLimNSubstrate <- CsynB <-	PhiB <- synB <- rG <- NA_real_
	computeAlphaFluxes <- function(alpha){
		ER <<- alpha * ETot
		EL <<- (1-alpha) * ETot
		#
		tvrER <<- alpha * synE
		tvrEL <<- (1-alpha) * synE
		#    
		limER <<- ER / (parms$kmR + ER)
		limEL <<- EL / (parms$kmL + EL)
		decR <<- decRp * limER
		decL <<- decLp * limEL
		#
		tvrERecycling <<- parms$kNB*(tvrER+tvrEL)
		uC <<- decC <<- decR + decL + tvrERecycling
		#
		#synE <<- parms$aEu * uC          # total enzyme production per uptake
		CsynBC <<- uC - synE/parms$eps - rM  # C required for biomass and growth respiration under C limitation
		#
		# Nitrogen balance
		decN <<- decR/cnR + decL/cnL + tvrERecycling/parms$cnE
		plantNUp <<- pmin(parms$plantNUp, decN/2)    # plants get at maximum half of the decomposed organic N
		PhiU <<- (1-parms$nu) * (decN-plantNUp)       # mineralization due to soil heterogeneity (Manzoni 08) 
		# immobilization flux
		immoPot <<- parms$iB * x["I"]
		uNSubstrate <<- (decN -plantNUp -PhiU)    # plant uptake also of organic N
		uNPot <<-  immoPot + uNSubstrate
		NsynBN <<- uNPot - synE/cnE
		CsynBN <<- (NsynBN*cnB)/parms$eps    # C required for biomass growth and growth respiration under N limitation
		#
		# compute C available for biomass, this time without taking into account immobilization flux
		NsynBNSubstrate <<- uNSubstrate - synE/cnE
		CsynBNSubstrate <<- (NsynBNSubstrate*cnB)/parms$eps    
		#
		isLimN <<- CsynBN <  CsynBC                      # N limited taking immobilization into account
		isLimNSubstrate <<-  CsynBNSubstrate < CsynBC    # N limited based on substrate uptake (without accounting for immobilization)
		CsynB  <<- if( isLimN ) CsynBN else CsynBC
		# average cN reqruired according to enzyme allocation accoording to C efficiencies
		#synE <- uC*parms$aE
		#uCGrowth <- uC - synE/parms$eps - rM    # C for both growth and growth respiration
		if( CsynB > 0){
			synB <<- parms$eps*CsynB
			rG <<- (1-parms$eps)*CsynB
		}else{
			synB <<- CsynB    # negative
			rG <<- 0
		}
		PhiB <<- uNSubstrate - (synE/parms$cnE + synB/parms$cnB) # imbalance mineralization/immobilization flux
		##value<< return the given argument
		alpha
	}
	alpha <- computeAlphaFluxes(alphaC)
	if( isLimN ){
		CsynBN_Clim = CsynBN
		CsynBC_Clim = CsynBC
		alpha <- computeAlphaFluxes(alphaN)
		if( !isLimN ){ 
			# when optimizing for C, system is in N-limitation
			# when optimizing for N, system is not in N-limitation
			# then comptute a balanced alpha
			#alpha <- computeAlphaFluxes((alphaN+alphaC)/2)
			alpha <- balanceAlphaSmallImmobilization(alphaC, alphaN, CsynBN_Clim, CsynBC_Clim
					, NsynBC=parms$eps*CsynBC/cnB, NsynBN)
			#c(alphaC=alphaC, alphaN=alphaN, alpha)
			#c(CsynBN_Clim/CsynBC_Clim, CsynBN/CsynBC)
		}
	}
	#c(alphaC, alphaN, alpha)
#recover()
	.tmp.f <- function(){
		alpha1 <- if( isTRUE(parms$useAlphaCUntilNLimited) || immoPot < uNSubstrate/100 ){ 
					balanceAlphaSmallImmobilization(alphaC, alphaN, CsynBN, CsynBC
							, NsynBC=parms$eps*CsynBC/cnB, NsynBN)
				} else {
					balanceAlphaLargeImmobilization(alphaC, alphaN, isLimN, isLimNSubstrate, immoAct=pmax(0,-PhiB), immoPot=immoPot)
				}
		computeAlphaFluxes(alpha1)
		# with new alpha new enzyme levels and changed PhiB, CsynBN, ...
		# do rudimentary fixpoint iteration
		alpha2 <- if( isTRUE(parms$useAlphaCUntilNLimited) || immoPot < uNSubstrate/100 ){ 
					balanceAlphaSmallImmobilization(alphaC, alphaN, CsynBN, CsynBC
							, NsynBC=parms$eps*CsynBC/cnB, NsynBN)
				} else {
					balanceAlphaLargeImmobilization(alphaC, alphaN, isLimN, isLimNSubstrate, immoAct=-PhiB, immoPot=immoPot)
				}
		alpha <- (alpha1+alpha2)/2 
		c(alphaC=alphaC, alphaN=alphaN, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha=alpha)
	}
    #
    # imbalance fluxes
    respSynE <- (1-parms$eps)/parms$eps * synE
    respO <- uC - (synE+respSynE+synB+rG+rM)
    MmImb <- uNPot - (synE/parms$cnE + synB/parms$cnB)      # imbalance taking into account potential immobilization
    #PhiB2 <- MmImb - immoPot # should equal PhiB
    PhiBU <- PhiB + PhiU
    respTvr <- (1-parms$epsTvr) * tvrB 
    PhiTvr <- respTvr/parms$cnB
    PhiBU <- PhiB + PhiU             # net microbial mineralization/immobilization flux when taking into account uptake mineralization
    PhiTotal <- PhiBU + PhiTvr       # total mineralization flux including microbial turnover
    #
    # Revenue strategy
    revLC <- decLp / (parms$kNL) / (parms$kmL + EL)
    revRC <- decRp / (parms$kNR) / (parms$kmR + ER)
    revLN <- revLC * cnE/cnL
    revRN <- revRC * cnE/cnR
    alphaC2 <- revRC / (revLC + revRC)  
    alphaN2 <- revRN / (revLN + revRN)
    #
    # tvr that feeds back to R pool, assume that N in SOM for resp (by epsTvr) is mineralized
    tvrC <- +parms$epsTvr*tvrB  +(1-parms$kNB)*(tvrER +tvrEL) 
    tvrN <- +parms$epsTvr*tvrB/parms$cnB  +(1-parms$kNB)*(tvrER +tvrEL)/parms$cnE 
    tvrExC <- 0     # fluxes leaving the system
    tvrExN <- 0
    # 
    leach <- parms$l*x["I"]
    #
    dB <- synB - tvrB
    dL <- -decL +parms$iL 
    dLN <- -decL/cnL  +parms$iL/parms$cnIL 
    dR <- -decR +parms$iR +tvrC 
    dRN <- -decR/cnR +parms$iR/parms$cnIR +tvrN 
    #dI <- +parms$iI +MmB +PhiTvr -(parms$kIP+parms$l)*x["I"]
    dI <- +parms$iI -parms$kIP -leach +PhiB +PhiU +PhiTvr        # plant uptake as absolute parameter
    #if( dI > 0.01 ) recover()    
    #
    if( isTRUE(parms$isFixedS) ){
        # scenario of fixed substrate
        dR <- dL <- dRN <- dLN <- dI <- 0
    }else if( isTRUE(parms$isTvrNil) ){
        # scenario of enzymes and biomass not feeding back to R
        dR <- +parms$iR -decR
        dRN <- +parms$iR/parms$cnIR -decR/cnR
        tvrExC <- tvrC
        tvrExN <- tvrN
    }
    respB <- respSynE + rG + rM + respO 
    resp <- respB + respTvr
    #
    resDeriv <-structure(as.numeric(c( dB, dR, dRN, dL, dLN, dI)),names=c("dB","dR","dRN","dL","dLN","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    sqrEps <- sqrt(.Machine$double.eps)
    # parms$iL - (decL + dL)
    # parms$iR + +tvrC -(decR + dR)
    # 
    if( diff( unlist(c(uC=uC, usage=respB+synB+synE )))^2 > sqrEps )  stop("biomass mass balance C error")
    #if( diff( unlist(c(uN=uNPot, usage=synE/parms$cnE + synB/parms$cnB + MmImb )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    if( diff( unlist(c(uN=uNSubstrate, usage=synE/parms$cnE + synB/parms$cnB + PhiB )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    if( !isTRUE(parms$isFixedS) ){    
        if( diff(unlist(c(dB+dR+dL+tvrExC+resp,    parms$iR+parms$iL )))^2 > sqrEps )  stop("mass balance C error")
        #if( diff(unlist(c( (dB)/parms$cnB+(dER+dEL)/parms$cnE+dRN+dLN+dI+tvrExN,    parms$iR/parms$cnIR +parms$iL/parms$cnIL -plantNUp +parms$iI -(parms$kIP+parms$l)*x["I"])))^2 > .Machine$double.eps )  stop("mass balance N error")
        if( diff(unlist(c( dB/parms$cnB +dRN+dLN+dI+tvrExN,    parms$iR/parms$cnIR +parms$iL/parms$cnIL -plantNUp +parms$iI -parms$kIP -parms$l*x["I"])))^2 > .Machine$double.eps )  stop("mass balance dN error")
    }
    #
    if( isTRUE(parms$isFixedR) ){ resDeriv["dR"] <- resDeriv["dRN"] <-  0   }        # for keeping R constant
    if( isTRUE(parms$isFixedL) ){ resDeriv["dL"] <- resDeriv["dLN"] <-  0   }        # for keeping L constant
    if( isTRUE(parms$isFixedI) ){ resDeriv["dI"] <-  0   }        # for keeping I constant
    #
    if( isTRUE(parms$isRecover) ) recover()    
    list( resDeriv, c(respO=as.numeric(respO)
		, ER=as.numeric(ER), EL=as.numeric(EL)
        , PhiB=as.numeric(PhiB), PhiU=as.numeric(PhiU), PhiTvr=as.numeric(PhiTvr) #, MmB=as.numeric(MmB)
        , PhiBU=as.numeric(PhiBU), PhiTotal=as.numeric(PhiTotal)
        , immoPot=as.numeric(immoPot), MmImb=as.numeric(MmImb)  
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
        , dR = as.numeric(dR), dL = as.numeric(dL), dB=as.numeric(dB), dI=as.numeric(dI)
    #    wCLim = (CsynBN/CsynBC)^delta
    #    wNLim = (parms$eps*CsynBC/cnB / NsynBN)^delta
        , uC=as.numeric(uC), synB=as.numeric(synB)
        , decC=as.numeric(decC), decN=as.numeric(decN)
))
}

.tmp.checkSteadyState <- function(){
	alphaC - decRp/(decRp + decLp*(parms$kN*parms$km+alphaC*parms$aE*x["B"])/(parms$kN*parms$km+(1-alphaC)*parms$aE*x["B"])) # check if really is a solution
	
	alphaC*synE - tvrER
	(1-alphaC)*synE - tvrEL
	revRC/(revLC+revRC) - alphaC
	decRp/(parms$kN*(parms$km+ER))
	decRp/(parms$kN*parms$km + alphaC*parms$aE*x["B"])
	
	parms$kL*decLp / (parms$kN*parms$km + (1-alphaC)*parms$aE*x["B"])
	parms$kR*decRp / (parms$kN*parms$km + (alphaC)*parms$aE*x["B"])
}



plotResSeam2 <- function(res, legendPos="topleft"
        , cls = c("B","ER","EL","respO","PhiTotal")
        , xlab="time (yr)"
        , ylab="gC/m2 , gN/m2 , % ,/yr"
        , subsetCondition=TRUE
        , ...
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
    res$I100 <- res$I  * 100 
    res$mPhiB10 <- -res$PhiB  * 10 
    res$mPhiBU10 <- -res$PhiBU  * 10 
    res$mPhiTotal10 <- -res$PhiTotal  * 10 
    #res$mMmB10 <- -res$MmB  * 10 
    #res$MmB100 <- res$MmB  * 100 
    #res$MmB1000 <- res$MmB  * 1000 
    res$cnR10 <- res$cnR  * 10
    res$dR <- c(NA, diff(res$R)/diff(res$time))
    res$dS <- c(NA, diff(res$R+res$L)/diff(res$time))
    #cls <- c("B100","ER_10","EL_10","respO","PhiTotal")
    #cls <- c("B100","ER_10","EL_10","respO","PhiTotal","eff")
    #cls <- c("B1000","ER_10","EL_10","PhiTotal")
    #cls <- c("B","E","respO","PhiTotal","R","L")
    matplot( res[subsetCondition,1], res[subsetCondition,cls], type="l", lty=1:20, col=1:20, xlab=xlab, ylab="", ...)
    mtext( ylab, 2, line=2.5, las=0)
    legend(legendPos, inset=c(0.01,0.01), legend=cls, lty=1:20, col=1:20)
}

