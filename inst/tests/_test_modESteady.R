#require(testthat)
context("modESteady")

test_that("modMeta setup correctly",{
            mm <- modMeta.ModESteady()
            expect_true( mm$nAux > 0)
            expect_equivalent(c("cf", "cs", "nf", "ns"), mm$colNames )
            expect_true(all(c("F", "S","B") %in% mm$rowNames) )
            namesMm <- c("matrixTemplate",   "elementNames",     "nCol"            
                    ,"nRow",            "colNames",         "csts"             
                    ,"rowNames",         "nAux",             "auxOutputTemplate"
                    ,"auxOutputNames",   "auxGroups", "iRUnits")
            expect_true( all(namesMm %in% names(mm)))
            expect_equivalent( mm$colNames, names(mm$iRUnits))
            expect_true( all(mm$iRUnits==1) )
        })


.setUp <- function(){
    .setUpDf <- within( list(),{
                modMeta <- modMeta.ModESteady();
                xc12=c(F=4,G=2,A=0.2,D=0.5,S=4)
                iR=c(s=1,a=2)
                modMeta$iRUnits["a"] = 0.1;
                (x0 <- initState.ModESteady(xc12=xc12, cn=0, iR=iR, modMeta=modMeta))
                #x0 <- modMeta$matrixTemplate; x0[] = 1:9
                
                xc12_noG=c(F=4,G=0,A=0.2,D=0.5,S=4)
                (x0_noG <- initState.ModESteady(xc12=xc12_noG, cn=0, iR=iR, modMeta=modMeta))
                
                xc12_noGDecay=c(F=0.004,G=0,A=2,D=0.5,S=4)
                (x0_noGDecay <- initState.ModESteady(xc12=xc12_noGDecay, cn=0, iR=iR, modMeta=modMeta))
                
                times <- seq(0,60)
                x0v <- structure(as.vector(x0), names=modMeta$elementNames)
            })
    attach(.setUpDf)
}

modMeta <- modMeta.ModESteady();
xc12=c(F=4,B=0.2,S=40)
parms <- within( list(),{
            a1=0.2
            a2=0.4
            betaE1 = 7
            betaE2 = 7
            betaB = 14
            mB = xc12["B"]/2            # half saturation for microbial biomass limitation
            kF=80/365.25  #; mFF=xc12["F"]		# Monod uptake of fresh material
            epsF=0.3; epsS=0.6			# maximum microbial efficiencies
            kS=0.1/365.25 #; mAS=xc12["A"]*10; eAS=1	# depolymerization with different equations
            sB=10/365.25						# maintenance respiration
            tvr=10/365.25; mmTvr=xc12["B"]*5; #tvrA=xc12["A"]*1.5; tvrP=1	#turnover by predation
            epsP=0.5; tvrS=0.1			# partitioning of turnover (includes predation, freezing, ...)
        })
p <- parms; p$modMeta <- modMeta;
cn = c(F=40,S=p$betaB, B=p$betaB)
amend = 3    
(x0 <- initState.ModESteady(xc12=xc12, cn=cn, modMeta=modMeta, amend=amend))


test_that("initialized correctly",{
    cn = c(F=40,S=p$betaB, B=p$betaB)
    amend = 3    
    (x1 <- initState.ModESteady(xc12=xc12, cn=cn, modMeta=modMeta, amend=amend))
    checkEquals( dimnames(modMeta$matrixTemplate), dimnames(x1) )
    expect_true(all(names(xc12) %in% modMeta$rowNames) )
    xc12 <- xc12[  modMeta$rowNames ]
    checkEqualsNumeric( xc12, x1[,"cs"] )
    checkEqualsNumeric( c(amend, 0,0), x1[,"ca"] )
    checkEqualsNumeric( xc12/cn, x1[,"ns"] )
    checkEqualsNumeric( c(amend/cn["F"], 0,0), x1[,"na"] )
})

checkOutputConsistency <- function(modMeta, aux){
    #aux <- res[[2]]
    checkEqualsNumeric( aux["respT"], sum(aux[c("respBiosynT","respOverflowT","respMaintT","respTvrT")]) )
    checkEqualsNumeric( aux["respT"], sum(aux[c("respFT","respGT","respAT","respDT","respTvrT")]) )
    checkEqualsNumeric( aux["respT"], sum( aux[paste("resp",modMeta$csts$c,sep="_")]*modMeta$iRUnits[modMeta$csts$c] ) )
}

checkMassBalance <- function(
        x0		##<< initial state
        ,out	##<< output of lsoda
        ,times	##<< times at which rates have been calculated
        ,iSteps	##<< integer vector of rows in out for which to check mass balance
        ,modMeta
){
    csums <- colSums(x0[,modMeta$csts$cis])
    respNames <- paste("resp",modMeta$csts$cis,sep="_")
    #dxNames <- paste("dxSums",modMeta$csts$cis,sep="_")
    resp <- apply( out[,respNames], 2, function(col){ twCalcObsScale(col)$oneSided} )
    #dx <- apply( out[,dxNames], 2, function(col){ twCalcObsScale(col)$oneSided} )
    iSteps <- iSteps[ iSteps!=0 & iSteps!=nrow(out) ]
    iSteps0 <- c(1,iSteps)
    stateMat <- twStateMatODERes(out, iSteps0, modMeta)
    checkEqualsNumeric( x0, stateMat[,,1] )
    #i=2
    #options(error=dump.frames)
    res <- sapply( seq_along(iSteps0), function(i){
                csumsState <- colSums(stateMat[,modMeta$csts$cis,i])
                csumsResp <- apply( resp, 2, twIntegrateTS, x=times, xr=times[ iSteps0[i] ] )
                #csumsDx <- apply( dx, 2, twIntegrateTS, x=times, xr=times[ iSteps0[i] ] )
                csumsi <- csumsState + csumsResp
                #checkEqualsNumeric(csums, csumsi, tolerance=sum(csums)*0.02)
                #checkEqualsNumeric(csums, csumsi, tolerance=sum(csums)*0.02)
                #c( csumsi, csumsResp, csumsDx )
                c( csumsi, csumsResp )
            })
    tmp.diff <- (res[1:2,] - csums)/csums
    expect_true( all(abs(tmp.diff)<1e-4) )	# relative error smaller than 1e-4
    #res
}


test_that("derivative is evaluated",{
    (x0 <- initState.ModESteady(xc12=xc12, cn=cn, modMeta=modMeta, amend=0))
    #trace(deriv.ModESteady, recover)
    res <- deriv.ModESteady(0,x0,p)
    checkEquals( dim(modMeta$matrixTemplate), dim(res[[1]]) )
    checkEquals( modMeta$auxOutputNames, names(res[[2]]) )
    #mtrace(checkOutputConsistency)
    checkOutputConsistency( modMeta, res[[2]] )
})

test.ode <- function(){
    #library(debug); mtrace(internal_test.ode)
    internal_test.ode(x0)
}

test.ode_noG <- function(){
    #library(debug); mtrace(internal_test.ode)
    internal_test.ode(x0_noG)
}

internal_test.ode <- function(x0){
    #mtrace(deriv.ModESteady)
    out <-  lsoda(x0, times, deriv.ModESteady, p )	
    
    colnames(out)
    y <- paste("csums",modMeta$rowNames,sep="_")
    matplot( out[,"time"], out[,y], type="l" )
    legend("topleft", y, col=1:6, lty=1:5, inset=0.02 )
    
    y2 <- c("respT", "respBiosynT", "respMaintT", "respTvrT")
    matplot( out[,"time"], out[,y2], type="l" )
    legend("topleft", y2, col=1:6, lty=1:5, inset=0.02 )
    
    auxArr <- (out[,-(1:(length(x0v)+1))])
    checkOutputConsistency(modMeta, auxArr[1])
    checkOutputConsistency(modMeta, auxArr[20])
    checkOutputConsistency(modMeta, auxArr[40])
    checkOutputConsistency(modMeta, auxArr[61])
    
    int <- floor(nrow(out) / 10.0)
    #mtrace(checkMassBalance)
    checkMassBalance( x0, out, times=times, iSteps=nrow(out)-(9:0)*int, modMeta=modMeta  )
}

test.derivC_decay <- function(){
    #mtrace(deriv.ModESteady)
    outR1 <- deriv.ModESteady(0,x0_noGDecay,p)	
    checkOutputConsistency(modMeta, outR1[[2]])
    expect_true( outR1[[2]]["actAT"] < 0)
    #outR1[[2]][ modMeta$auxGroups$tmp ]
    (out1 <- DLLfuncTest( x0_noGDecay, dllname = "twMDIHamer", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = p, times = 1,	nout = modMeta$nAux, outnames = modMeta$auxOutputNames)); 
    #tmp <- (out1[[2]] - outR1[[2]]); tmp[ abs(tmp) > .Machine$double.eps ]
    checkEqualsNumeric( outR1[[1]], out1[[1]] );
    #out1[[1]] - outR1[[1]] 
    checkEqualsNumeric( outR1[[2]], out1[[2]] );
    #with( as.list(out1[[2]]),c(decD_s,decD_a,tmp8 )) 
    #with( as.list(outR1[[2]]),c(kAppF_, (mFF_ + csums_F) ,csums_A / (mFF_ + csums_F), csums_A ))
}

checkDerivCParms <- function(p,x0){
    outR1 <- deriv.ModESteady(0,x0,p)	
    checkOutputConsistency(modMeta, outR1[[2]])
    #outR1[[2]][ modMeta$auxGroups$tmp ]
    
    #dynFilename <- file.path("src",paste("ModESteady", .Platform$dynlib.ext, sep = ""))
    #dynFilename <- file.path("src",paste("twMDIHamer", .Platform$dynlib.ext, sep = ""))
    #dyn.load(dynFilename); (out <- DLLfuncTest( x0v, dllname = "twMDIHamer", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = p, times = 1,	nout = modMeta$nAux, outnames = modMeta$auxOutputNames)); dyn.unload(dynFilename);
    (out1 <- DLLfuncTest( x0, dllname = "twMDIHamer", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = p, times = 1,	nout = modMeta$nAux, outnames = modMeta$auxOutputNames)); 
    #tmp <- (out1[[2]] - outR1[[2]]); tmp[ abs(tmp) > .Machine$double.eps ]
    checkEqualsNumeric( outR1[[1]], out1[[1]] );
    #out1[[1]] - outR1[[1]] 
    checkEqualsNumeric( outR1[[2]], out1[[2]] );
    #with( as.list(out1[[2]]),c(decD_s,decD_a,tmp8 )) 
    #with( as.list(outR1[[2]]),c(kAppF_, (mFF_ + csums_F) ,csums_A / (mFF_ + csums_F), csums_A ))
    
    
    #outR <-  lsoda(x0v, times, deriv.ModESteady, p )	
    outR <-  lsoda(x0, times, deriv.ModESteady, p )	
    auxArr <- (outR[,-(1:(length(x0)+1))])
    checkEquals( outR1[[2]], auxArr[1,])
    checkOutputConsistency(modMeta, auxArr[1,])
    checkOutputConsistency(modMeta, auxArr[20,])
    checkOutputConsistency(modMeta, auxArr[40,])
    checkOutputConsistency(modMeta, auxArr[61,])
    
    #checkEquals( modMeta.ModESteady(), modMeta ) FALSE adjusted modMeta$iRUnits["a"]
    #outR2 <- solve.ModESteady( x0, times, parms = p, modMeta=modMeta, useRImpl=TRUE);
    #checkEquals(outR, outR2)
    
    #dyn.load(dynFilename); out <- lsoda( x0v, times, dllname = "ModESteady", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = p, nout = modMeta$nAux, outnames = modMeta$auxOutputNames); dyn.unload(dynFilename);
    #out <- lsoda( x0v, times, dllname = "twMDIHamer", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = p, nout = modMeta$nAux, outnames = modMeta$auxOutputNames);
    #out <- solve.ModESteady( x0v, times, parms = p, useRImpl=FALSE);
    out <- solve.ModESteady( x0, times, parms = p, modMeta=modMeta, useRImpl=FALSE);
    # gives deviations in attributes checkEquals(outR, out)
    tmp <- outR-out
    .tmp.f <- function(){
        sort(abs(tmp[2,]),decreasing=TRUE)
        sort(abs(tmp[61,]),decreasing=TRUE)
    }
    tmpi <- which( tmp>.Machine$double.eps, arr.ind=TRUE)
    checkEqualsNumeric(rep(0,length(outR)), tmp )
    #str(attributes(outR))
    #str(attributes(out))
    
}

test.derivc <- function(){
    #library(debug); mtrace(checkDerivcParms)
    checkDerivCParms(p,x0)
}

test.derivc_noG <- function(){
    #library(debug); mtrace(checkDerivcParms)
    checkDerivCParms(p,x0_noG)
}

test.fSDec <- function(){
    ps <- p;
    idFGSDDec <- "fontaine"
    #mtrace(checkDerivcParms)
    for( idFGSDDec in names(decompFList) ){
        print(idFGSDDec)
        ps$fSDec <- idFGSDDec
        ps$f.sDec <- decompFList[[idFGSDDec]]
        checkDerivCParms(ps,x0)
    }
}

profile.ModESteady <- function(){
    #.setUp()
    runRmodN <- function (nRepeat,...){
        for( i in 1:nRepeat ){ outR <-  lsoda(x0v, times, deriv.ModESteady, p ) }
    }
    runCmodN <- function (nRepeat,...){
        #for( i in 1:nRepeat ){ out <- lsoda( x0v, times, dllname = "ModESteady", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = p, nout = modMeta$nAux, outnames = modMeta$auxOutputNames); }
        for( i in 1:nRepeat ){ out <- lsoda( x0v, times, dllname = "twMDIHamer", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = p, nout = modMeta$nAux, outnames = modMeta$auxOutputNames); }
    }
    n=10
    tR <- system.time( runRmodN(n) ) #2.6
    #dyn.load(dynFilename); tC <- system.time( runCmodN(n*100) ); dyn.unload(dynFilename); tC 
    tC <- system.time( runCmodN(n*100) ); #1.6 
    #speedup of about 160
    
    Rprof()
    tR <- system.time( runRmodN(n) ) #1.51
    Rprof(NULL)
    #?summaryRprof
    summaryRprof()
    
    Rprof()
    dyn.load(dynFilename); tC <- system.time( runCmodN(n*100) ); dyn.unload(dynFilename); tC 
    Rprof(NULL)
    #?summaryRprof
    summaryRprof()
    #seems that indeed the sum function over u is speeded up with C implementation
}

test.modVariantFS <- function(){
    # testing if the derivative with constrained parameters reprocudes FS model
    modMeta <- modMeta.ModESteady()
    #trace(parmsModelVariant.hamer, recover)
    if( 0==length(HamerParameterPriors$parms0$mmTvr) ) recover()
    #HamerParameterPriors$parms0$mmTvr
    parmsFS <- parmsModelVariant.hamer(c("tvr","obsBias"))
    xc12FS=c(F=4,G=0,A=0.2,D=0.5,S=4)
    iRFS=c(s=1,a=0)
    .amend=2.0
    (x0FS <- initState.ModESteady(xc12=xc12FS, cn=0, iR=iRFS, amend=.amend, usePrefSubst=FALSE))
    checkEqualsNumeric( xc12FS[modMeta$rowNames], x0FS[,"s"] )
    checkEqualsNumeric( .amend, x0FS["F","a"] )
    #mtrace(deriv.ModESteady)
    outR1<-NULL; outR1 <- deriv.ModESteady(0,x0FS,parmsFS$parms0)	
    #mtrace(checkOutputConsistency)
    checkOutputConsistency(modMeta, outR1[[2]])
    #outR1[[2]][ modMeta$auxGroups$tmp ]
    
    (out1 <- DLLfuncTest( x0FS, dllname = "twMDIHamer", func = "deriv_ModESteady",	initfunc = "init_ModESteady", parms = parmsFS$parms0, times = 1,	nout = modMeta$nAux, outnames = modMeta$auxOutputNames)); 
    #tmp <- (out1[[2]] - outR1[[2]]); tmp[ abs(tmp) > .Machine$double.eps ]
    checkEqualsNumeric( outR1[[1]], out1[[1]] );
    checkEqualsNumeric( outR1[[2]], out1[[2]] );
    
    outR <-  lsoda(x0FS, times, deriv.ModESteady, parmsFS$parms0 )	
    auxArr <- (outR[,-(1:(length(x0)+1))])
    checkEquals( outR1[[2]], auxArr[1,])
    checkOutputConsistency(modMeta, auxArr[1,])
    checkOutputConsistency(modMeta, auxArr[20,])
    checkOutputConsistency(modMeta, auxArr[40,])
    checkOutputConsistency(modMeta, auxArr[61,])
    
    out <- solve.ModESteady( x0FS, times, parms = parmsFS$parms0, modMeta=modMeta, useRImpl=FALSE);
    # gives deviations in attributes checkEquals(outR, out)
    bo <- outR[,"csums_A"] > 1e-5
    tmp <- (outR-out)
    .tmp.f <- function(){
        sort(abs(tmp[2,]),decreasing=TRUE)
        sort(abs(tmp[61,]),decreasing=TRUE)
    }
    tmpi <- which( tmp>.Machine$double.eps, arr.ind=TRUE)
    checkEqualsNumeric(rep(0,length(tmp[bo,])), as.vector(tmp[bo,]) )
    
    #check for same results as FS model
    xc12FSFS=c(F=4,A=0.2,S=4)
    (x0FSFS <- initState.SoilMod_FS(xc12=xc12FS, cn=0, iR=iRFS, amend=.amend))
    outFS <- solve.SoilMod_FS( x0FSFS, times, parms = HamerParameterPriors$parms0, modMeta=modMeta.SoilMod_FS(), useRImpl=FALSE);
    commonNames <- intersect(colnames(outR),colnames(outFS))	
    tmp <- outR[,commonNames]-out[,commonNames]
    .tmp.f <- function(){
        sort(abs(tmp[2,]),decreasing=TRUE)
        sort(abs(tmp[61,]),decreasing=TRUE)
    }
    tmpi <- which( tmp>.Machine$double.eps^0.5, arr.ind=TRUE)
    checkEqualsNumeric(rep(0,length(tmp)), tmp)
}

