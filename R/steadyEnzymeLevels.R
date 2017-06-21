computeAllocationPartitioning <- function(
	### compute allocation partitioning coefficient alpha assuming enzymes to be in quasi steady state
	dR		##<< potential decomposition flux from residue pool (gC or gN/time)
	,dL		##<< potential decomposition flux from labile pool
	,B		##<< microbial biomass carbon 
	,kM		##<< half-saturation constant in decomposition equation
	,kE		##<< enzyme turnover rate
	,aE		##<< proportion of microbial biomass allocated to enzyme production per time  
){
	c1 <- kM * kE
	aEB <- aE*B
	D <- 4*aEB^2*dL*dR + 8*aEB*c1*dL*dR + c1^2*dL^2 + 2*c1^2*dL*dR + c1^2*dR^2
	dLR <- dL - dR
    ans <- if( dLR == 0){
		0.5
	} else #if( dLR > 0 ){
		(-2*aEB*dR - c1*dL - c1*dR + sqrt(D))/(2*aEB*(dLR))
#	} else {
#		 -(2*aEB*dR + c1*dL + c1*dR + sqrt(D))/(2*aEB*(dLR))
#	}
	# ans - dR/(dR + dL*(kE*kM+ans*aE*B)/(kE*kM+(1-ans)*aE*B)) # check if really is a solution
	ans
}
attr(computeAllocationPartitioning,"ex") <- function(){
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
			,ER  = 2*parms0$km                  ##<< total enzyme pool
			,EL  = 4*parms0$km                   ##<< total enzyme pool
			,R = 1000                ##<< N rich substrate
			,RN = 1000/parms0$cnIR   ##<< N rich substrate N pool
			,L = 100                 ##<< N poor substrate
			,LN = 100/parms0$cnIL    ##<< N poor substrate N pool
			,I =  0.4                ##<< inorganic pool gN/m2
	)
	computeAllocationPartitioning( dR=parms$kR*x["R"], dL=parms$kL*x["L"], B=x["B"]
		,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	computeAllocationPartitioning( dR=parms$kR*x["R"], dL=parms$kL*x["L"]/500, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	computeAllocationPartitioning( dR=1, dL=1, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	computeAllocationPartitioning( dR=2, dL=1, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	computeAllocationPartitioning( dR=2, dL=0, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	computeAllocationPartitioning( dR=0, dL=1, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	computeAllocationPartitioning( dR=0, dL=0, B=x["B"]
			,kM = parms$km, kE=parms$kN, aE= parms$aE
	)
	
}