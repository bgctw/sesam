#  one-way:
# what is the steady state of E? E -> decomp -> uptake -> E prod -> E
# ran into cubic equation when tried to calculate steady state for given allocation to enzymes

modMeta.ModESteady <- function(
### Creating Meta-information for Basic Colimitation submodel.
){
    ##seealso<< \code{\link{twCreateModMeta}}, \code{\link{twSetModMetaAuxGroups}}
    modMeta <- twCreateModMeta(
            rowNames =  c('F','S','B')
            ,csts = list( cis=c('cs','ca'), nis = c('ns','na') )
    )
    auxGroups <- list(
            csums = paste("csums",modMeta$rowNames,sep="_")	#sum across all c compartments 
            ,decF = paste("decF",modMeta$colNames,sep="_") 
            ,decS = paste("decS",modMeta$colNames,sep="_")
            ,phi = "phi"
            ,tvrPred =  paste("tvrPred",modMeta$colNames,sep="_")
            ,resp = paste("resp",modMeta$csts$c,sep="_") 
            ,respT = "respT"	# total respiration over s and a
            # respiration by fueled process
            ,respBiosynT = "respBiosynT"
            ,respMaintT = "respMaintT"
            ,respTvrT = "respTvrT" 
            # respiration by compartments originating
            ,respFT = "respFT" 
            ,respST = "respST"
            ,tmp = paste("tmp",1:8,sep="")		#for debugging and getting information out of c-code
    )
    modMeta <- twSetModMetaAuxGroups(modMeta,auxGroups)
    #copy2clip(paste("enum AUX_OUTPUT_NAMES {",paste(c(modMeta.ModESteady()$auxOutputNames,"N_AUX"),collapse=","),"}; //generated in ModESteady.R",sep=""))	#to adjust soilmod_fs.c 
}


.depr.initState.ModESteady <- function(
        ### Creating initial state variables for Basic Colimitation submodel.
        xc12,	cn, 	iR 
        ,amend=0		##<< additional amendment of labelled carbon
        ,modMeta=modMeta.ModESteady()	##<< may pass pre-calulated modMeta for efficiency.
        ,usePrefSubst=TRUE	
        ### Wheter to use preferrential substrate utilization. 
        ### If false then amend goes to F instead of G pool. 
        ,... 
){
    ##seealso<< 
    x <- initStateModMeta(xc12,cn,iR,modMeta=modMeta)
    if( amend > 0 ) 
         x[ "F","f"] <- x[ "F","f"] + amend	
    x
    ### Numeric matrix (nPool, nIsotopes) of state variable mass.
}
#twUtestF("ModESteady",test="init")

initState.ModESteady <- function(
        ### Creating initial state variables for Basic Colimitation submodel.
        xc12,	cn 	 
        ,amend=0, cnAmend=cn["F"]		##<< additional amendment of labelled carbon
        ,modMeta=modMeta.ModESteady()	##<< may pass pre-calulated modMeta for efficiency.
        ### Wheter to use preferrential substrate utilization. 
        ### If false then amend goes to F instead of G pool. 
        ,... 
){
    ##seealso<< 
    x <- initStateModMeta(xc12,cn,iR=c(cs=1,ca=0,ns=1,na=0),modMeta=modMeta)
    if( amend > 0 ){ 
        x[ "F","ca"] <- x[ "F","ca"] + amend	
        x[ "F","na"] <- x[ "F","na"] + amend/cnAmend
    }
    x
    ### Numeric matrix (nPool, nIsotopes) of state variable mass.
}
#twUtestF("ModESteady",test="init")


solve.ModESteady <- function(
        ### solve the ODE of \code{\link{deriv.ModESteady}}
        x0		##<< numeric vector or matrix at t=0
        ,times	##<< times at which explicit estimates for y are desired. The first value in times must be the initial time.
        ,parms	##<< list of model parameters
        ,modMeta=modMeta.ModESteady()	##<< metaInformation from model. Pass for efficiency or when using different units. 
        ,useRImpl=FALSE	##<< flag indicating to use the R implementation instead of C implementation.
        ,useRk4=FALSE	##<< using the less precise but faster Runge-Kutta forth order function
        , useRImplMicCPart=FALSE
#,...			##<< only for R Implementation debugging purposes avoid ... to sort out errors 
){
    ##seealso<< \code{\link[deSolve]{lsoda}}
    # Tried to works with lsoda(..., atol=0). Else microbial biomass may become truely zero insteald of
    # small values, and derivative functions produces NaNs. But still gets 0
    parms$modMeta <- modMeta	#a way to pass it to the derivative function
    if( TRUE ){ #useRImpl ){ 
        #lsoda( x0, times, deriv.ModESteady, parms, atol = 0 )
        if( useRk4)
            rk4( x0, times, deriv.ModESteady, parms,  useRImplMicCPart=useRImplMicCPart)
        else
            lsoda( x0, times, deriv.ModESteady, parms,  useRImplMicCPart=useRImplMicCPart) 
    }else {
        #lsoda( x0, times, dllname = "twMDIHamer", func = "deriv_soilmod_fgsd3",	initfunc = "init_soilmod_fs", parms = parms, nout = modMeta$nAux, outnames = modMeta$auxOutputNames, atol = 0);
        if( useRk4)
            rk4( x0, times, dllname = "twMDIHamer", func = "deriv_soilmod_fgsd3",	initfunc = "init_soilmod_fgsd3", parms = parms, nout = modMeta$nAux, outnames = modMeta$auxOutputNames)
        else
            lsoda( x0, times, dllname = "twMDIHamer", func = "deriv_soilmod_fgsd3",	initfunc = "init_soilmod_fgsd3", parms = parms, nout = modMeta$nAux, outnames = modMeta$auxOutputNames)
    }
    ### \code{\link[deSolve]{lsoda}  
}

deriv.ModESteady <- function(
        ### Derivative function of Basic Colimitation model.
        t, x, p, useRImplMicCPart=FALSE
){
    #parms:
    # Meta information of this model passed with p
    mm <- p$modMeta
    # the derivative
    dx <- mm$matrixTemplate
    # auxiliary outputs 
    a <- mm$auxOutputTemplate
    
    # ode solvers provide x as a vector, attach attributes to have access as a mtrix.
    if( !is.matrix(x) ){ 
        dim(x) <- c(mm$nRow, mm$nCol)
        dimnames(x) <- dimnames(dx)
    }
    
    #carbon part
    csums <- sapply( mm$rowNames, function(rowName){  
                #sum( x[ rowName,mm$csts$cis] * unlist(mm$iRUnits[ mm$csts$cis ]) ) 	#take care that iRUnits holds all isotopes or parts
                sum( x[ rowName,mm$csts$cis] * mm$iRUnits[ mm$csts$cis ] ) 	#take care that iRUnits holds all isotopes or parts (in form c and not list for performance)
            }) #sum of all carbon fractions
    
    if( any(!is.finite(csums)) | (csums["B"] < sqrt(.Machine$double.eps)) ){
        # dead system: changes 0 as in matrixTemplate
        # dx is already 0
        a[ mm$auxGroups$csums ] = csums
        a[ mm$auxGroups$decF ] = 0
        a[ mm$auxGroups$decS ] = 0
        a[ mm$auxGroups$phi] = 0
        a[ mm$auxGroups$tvrPred ] = 0
        
        a[ mm$auxGroups$resp ] = 0
        a[ mm$auxGroups$respT ] = 0
        a[ mm$auxGroups$respBiosynT ] = 0
        a[ mm$auxGroups$respMaintT ] = 0
        a[ mm$auxGroups$respTvrT ] = 0
        
        a[ mm$auxGroups$respFT ] = 0
        a[ mm$auxGroups$respST ] = 0
    }else{
        #------ decomposition and uptake
        ##details<< 
        ## Assume assimilable to be in steady state: All material decomposed is taken up.
        ## In addition assume that Enzymes are in steady state: production of enzymes equals its turnover.
        #E1 <- p$tvrE1 
        limB <- csums["B"] / (p$mB + csums["B"])
        csumsUF <- csumsDecF <- p$kF * limB
        csumsUS <- csumsDecS <- p$kS * limB
        csumsU <- csumsDecF + csumsDecS
        decF <- uF <- if( csums["F"] > 0) csumsDecF * x["F",]/csums["F"] else 0*x["F",]
        decS <- uS <- if( csums["S"] > 0) csumsDecS * x["S",]/csums["S"] else 0*x["S",]
        u <- decF + decS
        #sum(u) - csumsU
        
        #------ internal partitioning
        ##details<<
        ## Maintenance respiration of active mics
        ## here not increasing with increasing A
        csumsRMaint = p$sB * csums["B"]
        phi <- (csumsU - csumsRMaint) / csumsU     # proportion to biosynthesis flux in uptake
        
        ##details<< 
        ## partitioning between biosynthesis and non-biosynthesis fluxes is proportional to 
        ## intrinsic SEU eps, i.e. the ratio of catabolic/anabolic fluxes for pure biosynthesis
        epsApp <- structure(phi * c(p$epsF, p$epsS), names=c("F","S")) # non-biosynthesis (maintenance) + biosynthesis respiration
        uAna <- u 
        uAna[mm$csts$cis] <- uF[mm$csts$cis]*epsApp["F"] + uS[mm$csts$cis]*epsApp["S"]
        csumsUAna <- sum( uAna[mm$csts$cis] * mm$iRUnits[ mm$csts$cis ])
        nsumsUAna <- sum( uAna[mm$csts$nis] * mm$iRUnits[ mm$csts$nis ])
        respInt <- u[mm$csts$cis] - uAna[mm$csts$cis]
        csumsRespInt <- sum( respInt * mm$iRUnits[ mm$csts$cis ])
        #(.respGr <- uF[mm$csts$cis]*phi*(1-p$epsF) + uS[mm$csts$cis]*phi*(1-p$epsS))
        #(.respMaint <- respInt - .respGr)
        #csumsRMaint - sum(.respMaint * mm$iRUnits[ mm$csts$cis ] )       # should equal zero
        #
        ##details<<
        ## Part of the biosynthesis is used for enzyme production.
        ## Different enzymes have different C/N ratios beta.
        ## The remainder is new biomass.
        csumsE1 <- p$a1 * csumsUAna
        csumsE2 <- p$a2 * csumsUAna
        nsumsE1 <- csumsE1 / p$betaE1
        nsumsE2 <- csumsE2 / p$betaE2
        uE1 <- p$a1 * uAna
        uE2 <- p$a2 * uAna
        uE1[mm$csts$nis] <- if( nsumsUAna > 0) nsumsE1 * uAna[mm$csts$nis]/nsumsUAna else 0*uAna[mm$csts$nis]
        uE2[mm$csts$nis] <- if( nsumsUAna > 0) nsumsE2 * uAna[mm$csts$nis]/nsumsUAna else 0*uAna[mm$csts$nis]
        #
        aB <- 1-(p$a1+p$a2)
        csumsB <- aB * csumsUAna
        nsumsB <- csumsB / p$betaB
        uB <- aB * uAna
        uB[mm$csts$nis] <- if( nsumsUAna > 0) nsumsB * uAna[mm$csts$nis]/nsumsUAna else 0*uAna[mm$csts$nis]
        #
        ##details<< 
        ## The nitrogen balance is closed by immobilization/mineralization flux.
        nsumsMin <- nsumsUAna - (nsumsE1 + nsumsE2 + nsumsB)
        nMin <- (uAna - (uE1 + uE2 + uB))[mm$csts$nis] 
        #
        #-------- external components
        # turnover by predation, rate increases with A
        #csumsPredA = p$tvr * csums["B"]^p$tvrP
        #csumsPredA = p$tvr * max(0,(csums["B"]-p$tvrA))^p$tvrP
        csumsPredA = p$tvr * csums["B"] * (csums["B"]/(csums["B"]+p$mmTvr))
        # isotopic signature of the microbes
        tvrPred <- csumsPredA/csums["B"] * x["B",]
        # distribute among respiration, and tvr to F and S
        tvrResp <- (1-p$epsP)*tvrPred
        tvrResp[mm$csts$nis] <- 0
        tvrPredF <- (1-p$tvrS)*(tvrPred - tvrResp)    
        tvrPredS <- (p$tvrS)*(tvrPred - tvrResp)
        # tvrPred - tvrResp - tvrPredF - tvrPredS   # all zero
        
        resp <- +respInt +tvrResp[mm$csts$cis]	 # total respiration
        dx["S",] <-  -decS +tvrPredS 
        dx["F",] <-  -decF +tvrPredF
        dx["B",] <-  +uB -tvrPred   
        
        #outputs
        a[ mm$auxGroups$csums ] = csums
        a[ mm$auxGroups$decF ] = decF[mm$colNames]
        a[ mm$auxGroups$decS ] = decS[mm$colNames]
        a[ mm$auxGroups$phi ] = phi
        a[ mm$auxGroups$tvrPred ] = tvrPred[mm$colNames]
        
        a[ mm$auxGroups$resp ] = resp[mm$csts$cis]
        a[ mm$auxGroups$respT ] = sum(resp[mm$csts$cis ] * mm$iRUnits[ mm$csts$cis ])
        a[ mm$auxGroups$respBiosynT ] = csumsRespInt - csumsRMaint
        a[ mm$auxGroups$respMaintT ] = csumsRMaint
        a[ mm$auxGroups$respTvrT ] = sum(tvrResp[mm$csts$cis]* mm$iRUnits[ mm$csts$cis ])
        
        a[ mm$tmp] = 0		
        
        #a[ mm$auxGroups$dxSums] = colSums(dx[,modMeta$csts$c])
        
        #a[ "kF_eff" ] = kF_eff["s"]
        if( any(!is.finite(dx))){	
            print(dx);
            dump.frames()
            recover()
            stop("produced non-finite values in derivative")	
        }
        if( any(!is.finite(a))){	
            print(dx);
            dump.frames()
            recover()
            stop("produced non-finite values in auxiliary outputs")	
        }
        # check mass balance
        tmp.balance = colSums(dx) + uE1 + uE2  #- colSums(tmp.input);
        tmp.balance[mm$csts$cis] <- tmp.balance[mm$csts$cis] + resp     
        tmp.balance[mm$csts$nis] <- tmp.balance[mm$csts$nis] + nMin
        if( any(abs(tmp.balance) > sqrt(.Machine$double.eps)) ){ 
            dump.frames()
            err.msg <- paste("mass balance violation: ",paste(tmp.balance,collapse=",")) 
            print(err.msg)
            recover()
            stop(msg)
        }
        
        #a[ mm$auxGroups$tmp[seq_along(mm$iRUnits)] ] <- mm$iRUnits
    }
    
    list(dx,a)
}




