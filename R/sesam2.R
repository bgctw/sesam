#library(deSolve)
# model corresponding to Seam2 with enzyme levels computed by quasi steady state 
# and neglecting mass fluxes by enzymes

# gC/m2 and gN/m2, /yr

derivSesam2 <- function(t,x,parms){
    x <- pmax(unlist(x),1e-16)      # no negative masses
    #
	# compute steady state enzyme levels for N and for C limitation
	decRp <- parms$kR * x["R"]
	decLp <- parms$kL * x["L"] 
	cnR <- x["R"]/x["RN"]
	cnL <- x["L"]/x["LN"]
	cnE <- parms$cnE #alphaC*parms$cnER + (1-alphaC)*parms$cnEL
	cnB <- parms$cnB
	alphaC <- computeSesam2AllocationPartitioning( dR=decRp, dL=decLp, B=x["B"]
			,kmkN = parms$kmkN, aE= parms$aE
	)
	alphaN <- computeSesam2AllocationPartitioning( dR=decRp/cnR, dL=decLp/cnL, B=x["B"]
			,kmkN = parms$kmkN, aE= parms$aE
	)
	#ETot <- parms$aE * x["B"] / parms$kN
	rM <- parms$m * x["B"]          # maintenance respiration
	tvrB <- parms$tau*x["B"]        # microbial turnover
	synE <- parms$aE * x["B"]       # total enzyme production per microbial biomass
	#respSynE <- (1-parms$eps)/parms$eps * synE	# growth respiration associated with enzyme production
	#
	# declare variables that will be computed/overidden in computeAlphaDependingFluxes
	# else operator '<<-' will override bindings in global environment
	decR <- decL <- uC <-   
			CsynBC <- decN <- plantNUp <- PhiU <- immoPot <- uNSubstrate <- uNPot <-
			NsynBN <- CsynBN <- isLimN <- 
			CsynB <- PhiB <- synB <- rG <- NA_real_
	computeAlphaDependingFluxes <- function(alpha){
		# compute fluxes that depend on alpha
		#tvrER <<- alpha * synE
		#tvrEL <<- (1-alpha) * synE
		#    
		decR <<- decRp * alpha*synE/(parms$kmkN + alpha*synE)
		decL <<- decLp * (1-alpha)*synE/(parms$kmkN + (1-alpha)*synE)
		#
		#tvrERecycling <<- parms$kNB*(tvrER+tvrEL)
		uC <<- decR + decL #+ tvrERecycling
		CsynBC <<- uC - rM #- synE/parms$eps # C required for biomass and growth respiration under C limitation
		#
		# Nitrogen balance
		decN <<- decR/cnR + decL/cnL #+ tvrERecycling/parms$cnE
		plantNUp <<- pmin(parms$plantNUp, decN/2)	# plants get at maximum half of the decomposed organic N
		PhiU <<- (1-parms$nu) * (decN-plantNUp)		# mineralization due to soil heterogeneity (Manzoni 08) 
		# immobilization flux
		immoPot <<- parms$iB * x["I"]
		uNSubstrate <<- (decN -plantNUp -PhiU)	# plant uptake also of organic N
		uNPot <<-  immoPot + uNSubstrate		# potential N uptake from both inorganic and organic
		NsynBN <<- uNPot #- synE/cnE			# N required for biomass growth (after invest into enzymes)
		CsynBN <<- (NsynBN*cnB)/parms$eps		# C required for biomass growth and associated growth respiration under N limitation
		#
		isLimN <<- CsynBN <  CsynBC				# wheter is microbial N limited (with taking immobilization into account)
		CsynB  <<- if( isLimN ) CsynBN else CsynBC
		if( CsynB > 0){
			synB <<- parms$eps*CsynB
			rG <<- (1-parms$eps)*CsynB
		}else{
			# with negative  biomass change, do not assign growth respiration
			synB <<- CsynB    
			rG <<- 0
		}
		PhiB <<- uNSubstrate - synB/parms$cnB #- synE/parms$cnE# imbalance mineralization/immobilization flux
		alpha	##value<< return the given argument
	}
	#
	alpha <- computeAlphaDependingFluxes(alphaC)		# first compute assuming C-limitation
	if( isLimN ){
		# if it turns out, that we are in N-limitation comptute again using alphaN
		CsynBN_Clim <- CsynBN; CsynBC_Clim = CsynBC # remember values for weights of C-limitation below
		alpha <- computeAlphaDependingFluxes(alphaN)
		if( !isLimN ){ 
			# co-limitation case:
			# when optimizing for C, system is in N-limitation
			# when optimizing for N, system is not in N-limitation
			# then comptute a balanced alpha
			alpha <- balanceAlphaBetweenCNLimitations(alphaC, alphaN, CsynBN_Clim, CsynBC_Clim
					, NsynBC=parms$eps*CsynBC/cnB, NsynBN)
			#c(alphaC=alphaC, alphaN=alphaN, alpha)
			#c(CsynBN_Clim/CsynBC_Clim, CsynBN/CsynBC)
		}
	}
    #
    # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
    respO <- uC - (synB+rG+rM) #-synE-respSynE
    respTvr <- (1-parms$epsTvr) * tvrB 
    PhiTvr <- respTvr/parms$cnB		# assuming same cn-Ratio of predators to equal cn ratio of microbes
    #
    # tvr that feeds back to R pool, assume that N in SOM for resp (by epsTvr) is mineralized
    tvrC <- +parms$epsTvr*tvrB  #+(1-parms$kNB)*(tvrER +tvrEL) 
    tvrN <- +parms$epsTvr*tvrB/parms$cnB  #+(1-parms$kNB)*(tvrER +tvrEL)/parms$cnE 
    tvrExC <- tvrExN <- 0     # fluxes leaving the system (will be set in scenario where turnover does not feed back)
    # 
    leach <- parms$l*x["I"]
    #
    dB <- synB - tvrB
    dL <- -decL +parms$iL 
    dLN <- -decL/cnL  +parms$iL/parms$cnIL 
    dR <- -decR +parms$iR +tvrC 
    dRN <- -decR/cnR +parms$iR/parms$cnIR +tvrN 
    dI <- +parms$iI -parms$kIP -leach +PhiU +PhiB +PhiTvr    # here plant uptake as absolute parameter
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
    #
    resDeriv <-structure(as.numeric(c( dB, dR, dRN, dL, dLN, dI)),names=c("dB","dR","dRN","dL","dLN","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    sqrEps <- sqrt(.Machine$double.eps)
    # parms$iL - (decL + dL)
    # parms$iR + +tvrC -(decR + dR)
    #
	# checking the mass balance of fluxes
	respB <- rG + rM + respO #+ respSynE 
	resp <- respB + respTvr
	#if( diff( unlist(c(uC=uC, usage=respB+synB+synE )))^2 > sqrEps )  stop("biomass mass balance C error")
	if( diff( unlist(c(uC=uC, usage=respB+synB )))^2 > sqrEps )  stop("biomass mass balance C error")
	#if( diff( unlist(c(uN=uNPot, usage=synE/parms$cnE + synB/parms$cnB + MmImb )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
    if( diff( unlist(c(uN=uNSubstrate, usage=synB/parms$cnB + PhiB )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
	#if( diff( unlist(c(uN=uNSubstrate, usage=synE/parms$cnE + synB/parms$cnB + PhiB )))^2 > .Machine$double.eps)  stop("biomass mass balance N error")
	if( !isTRUE(parms$isFixedS) ){    
        if( diff(unlist(c(dB+dR+dL+tvrExC+resp,    parms$iR+parms$iL )))^2 > sqrEps )  stop("mass balance C error")
        #if( diff(unlist(c( (dB)/parms$cnB+(dER+dEL)/parms$cnE+dRN+dLN+dI+tvrExN,    parms$iR/parms$cnIR +parms$iL/parms$cnIL -plantNUp +parms$iI -(parms$kIP+parms$l)*x["I"])))^2 > .Machine$double.eps )  stop("mass balance N error")
        if( diff(unlist(c( dB/parms$cnB +dRN+dLN+dI+tvrExN,    parms$iR/parms$cnIR +parms$iL/parms$cnIL -plantNUp +parms$iI -parms$kIP -parms$l*x["I"])))^2 > .Machine$double.eps )  stop("mass balance dN error")
    }
    #
	# allowing scenarios with holding some pools fixed
    if( isTRUE(parms$isFixedR) ){ resDeriv["dR"] <- resDeriv["dRN"] <-  0   }        # for keeping R constant
    if( isTRUE(parms$isFixedL) ){ resDeriv["dL"] <- resDeriv["dLN"] <-  0   }        # for keeping L constant
    if( isTRUE(parms$isFixedI) ){ resDeriv["dI"] <-  0   }        # for keeping I constant
	#
	# further computations just for output for tacking the system
	#ER <- alpha * parms$aE * x["B"] / parms$kN
	#EL <- (1-alpha) * parms$aE * x["B"] / parms$kN
	#limER <- ER / (parms$kmR + ER)
	#limEL <- EL / (parms$kmL + EL)
	revRC <- decRp / (parms$kmkN + alphaC*synE)
	revLC <- decLp / (parms$kmkN + (1-alphaC)*synE)
	revRN <- decRp/cnR / (parms$kmkN + alphaN*synE)
	revLN <- decLp/cnL / (parms$kmkN + alphaN*synE)
	PhiBU <- PhiB + PhiU             # net microbial mineralization/immobilization flux when taking into account uptake mineralization
	PhiTotal <- PhiBU + PhiTvr       # total mineralization flux including microbial turnover
	# c(alphaC, revRC/(revRC+revLC)); c(alphaN, revRN/(revRN+revLN)) # do not match in other limitation
	# compute C available for biomass, this time without taking into account immobilization flux
	NsynBNSubstrate <- uNSubstrate - synE/cnE
	CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
	isLimNSubstrate <-  CsynBNSubstrate < CsynBC    # N limited based on substrate uptake (without accounting for immobilization)
    #
    if( isTRUE(parms$isRecover) ) recover()    
    list( resDeriv, c(respO=as.numeric(respO)
		#, ER=as.numeric(ER), EL=as.numeric(EL)
        , PhiB=as.numeric(PhiB), PhiU=as.numeric(PhiU), PhiTvr=as.numeric(PhiTvr) #, MmB=as.numeric(MmB)
        , PhiBU=as.numeric(PhiBU), PhiTotal=as.numeric(PhiTotal)
        , immoPot=as.numeric(immoPot)  
        , alpha=as.numeric(alpha)
        , alphaC=as.numeric(alphaC), alphaN=as.numeric(alphaN)
        , cnR=as.numeric(cnR), cnL=as.numeric(cnL)
        #, limER=as.numeric(limER), limEL=as.numeric(limEL)
        , decR=as.numeric(decR), decL=as.numeric(decL)
        , resp=as.numeric(resp), respB=as.numeric(respB), respTvr=as.numeric(respTvr)
        , tvrB=as.numeric(tvrB)
        , revRC=as.numeric(revRC), revLC=as.numeric(revLC), revRN=as.numeric(revRN), revLN=as.numeric(revLN)
        , pCsyn = as.numeric(CsynBC / CsynBN), CsynReq=as.numeric(CsynBN), Csyn=as.numeric(CsynBC)
        , pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) ), NsynReq=as.numeric(CsynBC/cnB), Nsyn=as.numeric(NsynBN)
        , dR = as.numeric(dR), dL = as.numeric(dL), dB=as.numeric(dB), dI=as.numeric(dI)
        , uC=as.numeric(uC), synB=as.numeric(synB)
        , decN=as.numeric(decN)
))
}

computeSesam2AllocationPartitioning <- function(
		### compute allocation partitioning coefficient alpha assuming enzymes to be in quasi steady state
		dR		##<< potential decomposition flux from residue pool (gC or gN/time)
		,dL		##<< potential decomposition flux from labile pool
		,B		##<< microbial biomass carbon 
		,kmkN	##<< product of (half-saturation constant in decomposition equation) x (enzyme turnover rate)
		,aE		##<< proportion of microbial biomass allocated to enzyme production per time  
){
	aEB <- aE*B
	D <- 4*aEB^2*dL*dR + 8*aEB*kmkN*dL*dR + kmkN^2*dL^2 + 2*kmkN^2*dL*dR + kmkN^2*dR^2
	dLR <- dL - dR
	ans <- if( dLR == 0){
				0.5
			} else #if( dLR > 0 ){
				(-2*aEB*dR - kmkN*dL - kmkN*dR + sqrt(D))/(2*aEB*(dLR))
#	} else {
#		 -(2*aEB*dR + c1*dL + c1*dR + sqrt(D))/(2*aEB*(dLR))
#	}
	# ans - dR/(dR + dL*(kE*kM+ans*aE*B)/(kE*kM+(1-ans)*aE*B)) # check if really is a solution
	ans
}
attr(computeSesam2AllocationPartitioning,"ex") <- function(){
	parms <- list(
			cnB = 11
			,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
			#,cnIR = 4.5      ##<< between micr and enzyme signal
			,cnIR = 7      ##<< between micr and enzyme signal
			,cnIL = 30      ##<< N poor substrate, here no inorganic N supply, need to have low C/N for CO2 increase scenario
			,kN = 60       ##<< /yr enzyme turnover 60 times a year, each 6 days -> fast priming pulses
			,km = 0.05    ##<< enzyme half-saturation constant, in magnitude of enzymes, determined by kN
			,kR = 1/(10)        ##<< 1/(x years)       # to demonstrate changes on short time scale
			,kL = 5        ##<< 1/(x years)     # formerly 1 year
			,aE = 0.001*365   ##<< C biomass allocated to enzymes 1/day /microbial biomass
	)	
	x <- c( #aE = 0.001*365
			B = 17                     ##<< microbial biomass 
			,ER  = 2*parms$km                  ##<< total enzyme pool
			,EL  = 4*parms$km                   ##<< total enzyme pool
			,R = 1000                ##<< N rich substrate
			,RN = 1000/parms$cnIR   ##<< N rich substrate N pool
			,L = 100                 ##<< N poor substrate
			,LN = 100/parms$cnIL    ##<< N poor substrate N pool
			,I =  0.4                ##<< inorganic pool gN/m2
	)
	computeSesam2AllocationPartitioning( dR=parms$kR*x["R"], dL=parms$kL*x["L"], B=x["B"]
			,kmkN = parms$km*parms$kN, aE= parms$aE
	)
	computeSesam2AllocationPartitioning( dR=parms$kR*x["R"], dL=parms$kL*x["L"]/500, B=x["B"]
			,kmkN = parms$km*parms$kN, aE= parms$aE
	)
	computeSesam2AllocationPartitioning( dR=1, dL=1, B=x["B"]
			,kmkN = parms$km*parms$kN, aE= parms$aE
	)
	computeSesam2AllocationPartitioning( dR=2, dL=1, B=x["B"]
			,kmkN = parms$km*parms$kN, aE= parms$aE
	)
	computeSesam2AllocationPartitioning( dR=2, dL=0, B=x["B"]
			,kmkN = parms$km*parms$kN, aE= parms$aE
	)
	computeSesam2AllocationPartitioning( dR=0, dL=1, B=x["B"]
			,kmkN = parms$km*parms$kN, aE= parms$aE
	)
	computeSesam2AllocationPartitioning( dR=0, dL=0, B=x["B"]
			,kmkN = parms$km*parms$kN, aE= parms$aE
	)
}