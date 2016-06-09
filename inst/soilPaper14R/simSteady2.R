# using Seam2 version with nitrogen efficiency nu
isPaperBGC <- ("paperBGC" %in% commandArgs(trailingOnly = TRUE))
#isPaperBGC <- TRUE

# simulating CO2 increase and bare soil, based on modEezy5
baseFontSize <- 16  # presentations
# based on 
if( isPaperBGC){
    library(twDev)
    loadPkg()
    baseFontSize <- 9  # pubs
} 
library(ggplot2)
library(grid)   #unit
library(reshape)  # melt
library(RColorBrewer) # brewer.pal

# use the same colors in every graph - color-blind friendly
myColors <- brewer.pal(5,"Dark2")
.simpleCap <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2),  sep = "")
}
#names(myColors) <- .simpleCap(names(parmsScen))
colScale <- scale_colour_manual(name = "Allocation",values = myColors)


# gC/m2 and gN/m2, /yr
parms0 <- list(
        cnB = 7.16
        ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        #,cnIR = 4.5      ##<< between micr and enzyme signal
        ,cnIR = 7      ##<< between micr and enzyme signal
        ,cnIL = 30      ##<< N poor substrate, here no inorganic N supply, need to have low C/N for CO2 increase scenario
        #,cnIL = 40      ##<< N poor substrate       # near colimitation, here N limited (no overflow)
        #,cnIL = 50      ##<< N poor substrate       # near colimitation, here N limited (no overflow)
        #,kN= 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
        #,kN= 1/(1/12)  ##<< 1 month (Blagodatskaya 60 days priming)
        #,km = 0.05     ##<< enzyme half-saturation constant
        ,kN = 60       ##<< /yr enzyme turnover 60 times a year, each 6 days -> fast priming pulses
        ,km = 0.01    ##<< enzyme half-saturation constant, in magnitude of enzymes, determined by kN
        ,kNB = 0.8      ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R)
        #,kR = 1/(50)        ##<< 1/(x years) 
        ,kR = 1/(10)        ##<< 1/(x years)       # to demonstrate changes on short time scale
        #,kR = 1/(20)        ##<< 1/(x years)       # to demonstrate changes on short time scale
        #,kL = 1/(0.33)        ##<< 1/(x years)     # formerly 1 year
        ,kL = 5        ##<< 1/(x years)     # formerly 1 year
        ,aE = 0.001*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,m = 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = 1/60*365  ##<< biomass turnover rate (12 days)
        ,eps = 0.5      ##<< carbon use efficiency
        ,epsTvr = 0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
        ,iR = 0        ##<< input modelled explicitely
        #,iL = 300         # g/m2 input per year (half NPP)
        ,iL = 400         # g/m2 input per year (half NPP)
        #,plantNUp = 300/70*1/4  # plant N uptake balancing N inputs
        ,plantNUp = 0
        ,useFixedAlloc=FALSE    ##<< set to true to use Fixed enzyme allocation (alpha = 0.5)
        ,iP = 10.57 #0.0289652*365          ##<< plant uptake iP I
        ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I
        ,iI = 0     ##<< input of mineral N
        ,l = 0.96   #0.00262647*365       ##<< leaching rate of mineralN l I

        ,kIP = 10.57 #0.0289652*365          ##<< plant uptake iP I
        ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I
        ,iI = 0     ##<< input of mineral N
        ,l = 0   #0.00262647*365       ##<< leaching rate of mineralN l I
        ,nu = 1     # microbial N use efficiency
        ,useAlphaCUntilNLimited = TRUE      ##<< do not decrease investment into C enzmyes when NSubstrateLimtited, but only when N-Limited
)

parms0 <- within(parms0,{
            kmR <- kmL <- km
            epsR <- epsL <- eps
            cnER <- cnEL <- cnE 
            kNR <- kNL <- kN
        })
parms <- parms0

x0 <- x0Orig <- c( #aE = 0.001*365
        B = 17                     ##<< microbial biomass 
        ,ER  = 2*parms0$km                  ##<< total enzyme pool
        ,EL  = 4*parms0$km                   ##<< total enzyme pool
        ,R = 400                 ##<< N rich substrate
        ,RN = 400/parms0$cnIR    ##<< N rich substrate N pool
        ,L = 100                 ##<< N poor substrate
        ,LN = 100/parms0$cnIL    ##<< N poor substrate N pool
        ,I =  0                  ##<< inorganic pool gN/m2
)
x <- x0


parmsScen <- list(
        Revenue = parms0
        ,Fixed = within(parms0, {isAlphaFix <- TRUE})
        ,Match = within(parms0, {isAlphaMatch <- TRUE})
)
names(myColors) <- names(parmsScen)




simfCNGraph <- function(
    ### fixing Substrate availability and varying N amount in L
){
    parms0F <- within(parms0,{    isFixedS <- TRUE; kNB = 0   })
    parmsScenF <- list(
            Revenue = parms0F
            ,Fixed = within(parms0F, {isAlphaFix <- TRUE})
            ,Match = within(parms0F, {isAlphaMatch <- TRUE})
    )
    
    cnR <- 6.8
    cnL <- 30
    x0N <- c( #aE = 0.001*365
            B = 20                     ##<< microbial biomass 
            ,ER  = 1.5*parms0$km                  ##<< total enzyme pool
            ,EL  = 1.5*parms0$km                   ##<< total enzyme pool
            ,R = 400                 ##<< N rich substrate
            ,RN = 400/cnR    ##<< N rich substrate N pool
            ,L = 100                  ##<< N poor substrate
            ,LN = 100/cnL    ##<< N poor substrate N pool
            ,I =  0.1                    ##<< inorganic pool
    )
    
    cnLs <- seq( 18,62,by=1)
    cnLs <- seq( 18,40,by=2)
    #cnLs <- seq( 18,80,by=0.5)
    #cnLs <- seq( 18,160,by=5)
    cnL <- 23
    resLC <- lapply( cnLs, function(cnL){ 
                cat(cnL,", ")
                x0N["LN"] <- x0N["L"]/cnL
                times <- seq(0,10, length.out=101)
                #times <- seq(0,10000, length.out=101)
                #scen <- "Match"
                #scen <- "Revenue"
                resL <- lapply( names(parmsScenF), function(scen){
                #resL <- lapply( "Revenue", function(scen){
                            parmsInit <- parmsScenF[[scen]]
                            res <- res1 <- as.data.frame(lsoda( x0N, times, derivSeam2, parms=parmsInit))
                            xE <- unlist(tail(res,1))
                            #plotResSeam1(res, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
                            #plotResSeam1(res, "topright", cls = c("I"))
                            #trace(derivSeam2, recover) #untrace(derivSeam2)
                            #tmp <- derivSeam2(0, xE[1:length(x0)+1], parmsInit)
                            xE
                        })
                #lapply(resL, "[[", "alpha")
                resScen <- cbind( data.frame(scen=names(parmsScenF), do.call( rbind, resL )))
            })
    resC <- resCAll <- do.call( rbind, resLC )
    #res <- resLC[[7]] 
    resC <- subset(resCAll, cnL <= 160)
    resC$MmImbMon <- resC$MmImb/12
    resC$respOMon <- resC$respO/12

    #p8a <- ggplot(subset(resC, scen=="Revenue"), aes(x=cnL, y=isLimN)) + geom_line()
    #p8a
    resCRev <- subset(resC, scen=="Revenue")
    cnTER <- resCRev$cnL[ which.min( abs(1-resCRev$pNsyn) ) ]
    
    dsp <- melt(resC, id=c("cnL","scen"), measure.vars=c("alpha","B","MmImbMon","respOMon"), variable_name="Measure")
    #levels(dsp$Measure) <- c("Allocation_Ratio (alpha)","Biomass (gC/m^2)","Mineralization[Imb] (gN/m^2/yr)","OverflowRespiration (gC/m2/yr)")
    levels(dsp$Measure) <- c("Allocation~to~ R~(alpha)","Biomass~(gCm^{-2})","Min~Imb~(gNm^{-2}*month^{-1})","Overflow~(gCm^{-2}*month^{-1})")
    p8 <- ggplot( dsp, aes(x=cnL, y=value, col=scen), environment = environment() ) + geom_line(size=1) + 
            #facet_grid(Measure~., scales = "free", labeller= label_parsed ) + 
            facet_wrap(~Measure, scales = "free_y", labeller = label_parsed ) +
            scale_x_continuous('C/N ratio of Litter') + 
            geom_vline(aes(xintercept=cnTER), colour="#990000", linetype="dashed") +
            theme_bw(base_size=baseFontSize) +
            theme(axis.title.y = element_blank()) +
            theme()
    #twWin(6)
    p8 + colScale
    
    if (isPaperBGC){
        twWin(width=3.3, height=3.2, pointsize=9, pdf="soilPaper14/fig/VarNNoFeedback.pdf")
        print( p8 + colScale + theme(legend.position = c(0.14,0.475), legend.justification=c(0,1)) ) 
        dev.off()
    }

    dsp <- melt(resC, id=c("cnL","scen"), measure.vars=c("alpha","alphaC","alphaN","pCsyn", "pNsyn"), variable_name="Measure")
    #levels(dsp$Measure) <- c("Allocation_Ratio (alpha)","Biomass (gC/m^2)","Mineralization[Imb] (gN/m^2/yr)","OverflowRespiration (gC/m2/yr)")
    #levels(dsp$Measure) <- c("Allocation to E_SOM (alpha)","Biomass (gC/m2)","Mineralization Imb (gN/m2/month)","OverflowRespiration (gC/m2/month)")
    p8a <- ggplot( dsp, aes(x=cnL, y=value, col=scen), environment = environment() ) + geom_line(size=1) + 
            #facet_grid(Measure~., scales = "free", labeller= label_parsed ) + 
            facet_wrap(~Measure, scales = "free_y" ) +
            #facet_wrap(~Measure ) +
            scale_x_continuous('C/N ratio of Lit') + 
            geom_vline(aes(xintercept=cnTER), colour="#990000", linetype="dashed") +
            theme_bw(base_size=baseFontSize) +
            theme(axis.title.y = element_blank()) +
            theme()
    #twWin(6)
    p8a 
    
    
    dspC <- melt(subset(resC, scen=="Revenue"), id=c("cnL"), measure.vars=c("revRC","revLC"), variable_name="Pool"); dspC$Element = "C"
    levels(dspC$Pool) <- c("R","L")
    dspN <- melt(subset(resC, scen=="Revenue"), id=c("cnL"), measure.vars=c("revRN","revLN"), variable_name="Pool"); dspN$Element = "N"
    levels(dspN$Pool) <- c("R","L")
    dsp <- rbind(dspC, dspN)
    #levels(dsp$Measure) <- c("Allocation_Ratio (alpha)","Biomass (gC/m^2)","Mineralization[Imb] (gN/m^2/yr)","OverflowRespiration (gC/m2/yr)")
    #levels(dsp$Measure) <- c("Allocation ratio (alpha)","Biomass (gC/m2)","Mineralization Imb (gN/m2/month)","OverflowRespiration (gC/m2/month)")
    p8b <- ggplot( dsp, aes(x=cnL, y=value, col=Pool), environment = environment() ) + geom_line(size=1) + 
            #facet_grid(Measure~., scales = "free", labeller= label_parsed ) + 
            facet_wrap(~Element, scales = "free_y" ) +
            scale_x_continuous('C/N ratio of Lit') +
            geom_vline(aes(xintercept=cnTER), colour="#990000", linetype="dashed") +
            theme_bw(base_size=baseFontSize) +
            #theme(axis.title.y = element_blank()) +
            ylab("Revenue") +
            theme()
    #twWin(6)
    p8b 
}



simInitSteady <- function(
    ### inspect approaching a steady state (or breakdown of biomass)
){
    scen <- "Match"
    scen <- "Fixed"
    scen <- "Revenue"
    resAll <- lapply( c("Revenue","Fixed","Match"), function(scen){
                parmsInit <- parmsScen[[scen]]
                parmsInit$isFixedI <- TRUE
                times <- seq(0,150, length.out=301)
                #times <- seq(0,10000, length.out=101)
                #res <- res1 <- as.data.frame(lsoda( x0, times, derivSeam2, parms=within( parmsInit, isRecover <- TRUE)))
                res <- res1 <- as.data.frame(lsoda( x0, times, derivSeam2, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivSeam2, parms=within(parms0, useFixedAlloc<-TRUE) ))
                xE <- unlist(tail(res,1))
                #plotResSeam1(res, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
                #plotResSeam1(res, "topright", cls = c("ER","EL"))
                tail(res)
                #head(res[,c("limE1","limE2")])
                #tail(res[,c("limE1","limE2")])
                res$scen <- scen
                res
            })
    resScen <- do.call( rbind, resAll)
    
    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("R","L"),variable_name="Pool")
    names(dsp)[names(dsp)=="scen"] <- "Allocation"
    levels(dsp$Pool) <- c("R","L")
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area() 

    .tmp.f <- function(){
        p1 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
                theme_bw(base_size=baseFontSize) +
                ylim(c(100,3000)) +
                labs(x="Time (yr)", y="Carbon stock (gC/m2)", linetype="Substrate pool") +
                theme()                
        p1+ colScale
    }
   
    dsp$value[ dsp$value > 700] <- NA
    p1b <- ggplot( dsp, aes(x=time, y=value, col=Allocation)) + geom_line(size=1) +
            facet_wrap( ~ Pool,scales="free_y") + 
            theme_bw(base_size=baseFontSize) +
            labs(x="Time (yr)", y="Carbon stock (gC/m2)", linetype="Substrate pool") +
            theme()                
    p1b+ colScale
    
    
    if (isPaperBGC){
        twWin(width=3.3, height=2, pointsize=9, pdf="soilPaper14/fig/SimSteady.pdf")
        print(p1b + colScale + theme(legend.position = c(1.05,1.05), legend.justification=c(1,1))
                )
        dev.off()
    }
    
    
}

               
simCO2Increase <- function(
    ### Simulated increase of C-input by 20% during years 10-60
){
    #scen <- "Revenue"
    t1S <- 10
    t2I = 50
    t3S <- 50
    fInputInc = 1.2
    resAll <- lapply( c("Revenue","Fixed"), function(scen){
            parmsInit <- parmsScen[[scen]]
            parmsInit$isFixedI <- TRUE
            #parmsInit$cnIL <- 30
            x0Pr <- x0
            #x0Pr["I"] <- 0.8
            x0Pr["R"] <- 1100
            x0Pr["RN"] <- x0Pr["R"]*parmsInit$cnIR
            # spinup run
            times <- seq(0,500, length.out=101)
            #times <- seq(0,10000, length.out=101)
            res <- res1 <- as.data.frame(lsoda( x0Pr, times, derivSeam2, parms=parmsInit))
            #res <- res1f <- as.data.frame(lsoda( x0, times, derivSeam2, parms=within(parms0, useFixedAlloc<-TRUE) ))
            xE <- unlist(tail(res,1))
            plotResSeam1(res, "topright", cls = c("B10","respO","MmB","Rr","Lr","alpha100"))
            #plotResSeam1(res, "topright", cls = c("ER","EL"))
            #tmp <- derivSeam2(0, xE[1:length(x0)+1], within(parmsInit, isRecover <- TRUE) )
                
            
            # 10 yr steady state
            res <- res1S <- as.data.frame(lsoda( xE[1:length(x0)+1], 1:t1S, derivSeam2, parms=parmsInit))
            
            # 30 yr double C input    
            parmsC2 <- within(parmsInit, { # more C in litter, but not more N 
                        iL <- iL*fInputInc
                        cnIL <- cnIL*fInputInc
                        #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
                    }) 
            times <- seq(0,t2I, length.out=101)
            res <- res2I <- as.data.frame(lsoda( xE[1:length(x0)+1], times, derivSeam2, parms=parmsC2))
            res2I$time <- res2I$time +t1S 
            xE2 <- unlist(tail(res2I,1))
            plotResSeam1(res, "topright", cls = c("B10","respO","MmB","Rr","Lr","alpha100"))
            #plotResSeam1(res, "topright", cls = c("I"))
            #plotResSeam1(res, "topright", cls = c("R","L"))
            #trace(derivSeam2, recover)        #untrace(derivSeam2)
            #tmp <- derivSeam2(0, xE[1:length(x0)+1], within(parmsC2, isRecover <- TRUE) )
            #tmp <- derivSeam2(0, xE2[1:length(x0)+1], within(parmsC2, isRecover <- TRUE) )
            
            # 30 yr normal input
            times <- seq(0,t3S, length.out=101)
            res <- res3S <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivSeam2, parms=parmsInit))
            res3S$time <- res3S$time +t1S + t2I
            xE3 <- tail(res3S,1)
            #plotResSeam1(res, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
            
            resc <- rbind( res1S, res2I[-1,], res3S[-1,])
            resc$scen <- scen
            #plotResSeam1(resc, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
            resc
    })
    resScen <- do.call( rbind, resAll)

    resScen$RL <- resScen$R + resScen$L
    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("L","R","RL"),variable_name="Pool")
    levels(dsp$Pool) <- c("L","R","L+R")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    #names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area()

    .tmp.f <- function(){
        p2 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
                xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") + labs(linetype="Substrate pool") +
                theme_bw(base_size=baseFontSize) +
                #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
                theme()                
        p2 + colScale
    }

    p2f <- ggplot( dsp, aes(x=time, y=value, col=Allocation)) + geom_line(size=1) + 
            facet_grid(Pool ~ .,scales="free_y") + 
            xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") + labs(linetype="Substrate pool") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p2f + colScale
    
    if (isPaperBGC){
        twWin(width=3.3, height=3.3, pointsize=9, pdf="soilPaper14/fig/CO2Increase.pdf")
        print(p2f + colScale+ theme(legend.position = c(1.0,1.0), legend.justification=c(1,1)) +theme( plot.margin = unit( c(0,0,0,0) , "in" ) )
        )
        dev.off()
    }
    
    
    # graph of inputs
    dsCN <- data.frame( time=1:(t1S+t2I+t3S), iL=parmsScen[["Revenue"]]$iL )
    dsCN[ dsCN$time %in% (t1S+1):(t1S+t2I),"iL"] <- parmsScen[["Revenue"]]$iL * fInputInc

    p2b <- ggplot( dsCN, aes(x=time, y=iL)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Input") + ylim(c(0,500)) +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()
    p2b
    
    # imbalance fluxes
    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("respO","MmImb"),variable_name="Pool")
    levels(dsp$Pool) <- c("Overflow respiration","N Mineralization")
    #dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("resp","respO","Mm","MmImb","B","RN"),variable_name="Pool")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    #names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area()
    
    p2 <- ggplot( dsp, aes(x=time, y=value, col=Allocation)) + geom_line(size=1) +
            facet_grid(Pool ~ .,scales="free_y") + 
            xlab("Time (yr)")+ ylab("Imbalance flux (g(C or N)/m2/yr)") + labs(linetype="Ouput fluxes") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme(legend.position = c(0.95,1.0), legend.justification=c(1,1)) +
            theme()                
    p2 + colScale

    if (isPaperBGC){
        twWin(width=3.3, height=2.35, pointsize=9, pdf="soilPaper14/fig/CO2IncreaseImb.pdf")
        print(p2 + colScale)  
        dev.off()
    }
    
    
    .tmp.f <- function(){
        dsCN$graph <- "Litter C Input"
        dsp2 <- dsp
        dsp2$graph <- "Substrate Pools"
        dspM <- rbind( dsCN)
    }
}
    



simPriming <- function(
### Simulated decrease of input to 1/100 during years 10-210
){
    #scen <- "Fixed"
    #scen <- "Revenue"
    resAll <- lapply( c("Revenue","Fixed"), function(scen){
                parmsInit <- parmsScen[[scen]]
                parmsInit$isFixedI <- TRUE
                #parmsInit$kNL <- parmsInit$kNR <- 60; parmsInit$kmL <- parmsInit$kmR <- 0.005
                #parmsInit$kR <- 1/5 # tested whether this speeds up dynamics, but it does not 
                parmsInit$kL <- (365/0.1) # 10th day (unlimited) turnover time of L
                
                #parmsInit$kL = 1/(1/365)       # use a easily degradable substrate               
                # sim steady state
                times <- seq(0,500, length.out=20)
                #times <- seq(0,500, length.out=101)
                res <-  as.data.frame(lsoda( x0, times, derivSeam2, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivSeam2, parms=within(parms0, useFixedAlloc<-TRUE) ))
                #plotResSeam1(res, "topright", cls = c("B10","respO","L","R","alpha100"))
                xESteadyP <- xE <- tail(res,1)
                xE
                
                # X yr decreased C input    
                t2I = if( scen=="Revenue") 18 else 10   # revenue needs a bit longer, make sure about the same base level
                fInputInc = 1.2
                parmsC2 <- within(parmsInit, {  
                            #iL <- iL/10
                            iL <- iL/100
                            #cnIL <- cnIL*fInputInc
                            #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
                        }) 
                times <- c(seq(0,t2I, length.out=101))
                res <- as.data.frame(lsoda( unlist(xESteadyP[1:length(x0)+1]), times, derivSeam2, parms=parmsC2))
                res2I <- tail(res,10)   # last 10 years
                xE2 <- unlist(tail(res2I,1))
                #plotResSeam1(res, "topright", cls = c("B10","respO","MmB","Rr","Lr","alpha100"))
                #-----  some long-time fluctutions of break-out of biomass: emergent!
                #plotResSeam2(res2I, "topright", cls = c("B10","respO","MmB","Rr","Lr","alpha100"))
                #plotResSeam2(res, "topright", cls = c("L"))
                #trace(derivSeam2, recover)        #untrace(derivSeam2)
                #tmp <- derivSeam2(0, xE2[1:length(x0)+1], within( parmsC2  ,isRecover <-TRUE))
                
                # 5 yr control of continued decreased input 
                t3S <- 2
                #t3S <- 5
                times <- seq(0,t3S, length.out=501)
                res <- res3Sc <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivSeam2
                                , parms=parmsC2
                                #, parms=within( parmsC2  ,isRecover <-TRUE)
                        ))
                #plotResSeam2(res3Sc, "topright", cls = c("B10","respO","Rr","Lr","alpha100","EL_1000","ER_1000","Phi","Mm"))
                res3Sc$time <- res3Sc$time 
                xE3c <- tail(res3Sc,1)
                
                # 5 yr of continued decreased input but with initial pulse of L
                x0P <- xE2; x0P["L"] <- x0P["L"] + 50 
                res <- res3S <- res3Sp <- as.data.frame(lsoda( x0P[1:length(x0)+1], times, derivSeam2
                                        , parms=parmsC2
                                        #, parms=within( parmsC2  ,isRecover <-TRUE)
                ))
                #plotResSeam2(res3Sp, "topright", cls = c("B10","respO","Rr","Lr","alpha100","EL_1000","ER_1000"))
                #plotResSeam2(res3Sp, "topright", cls = c("B10","Lr","alpha100","EL_1000","ER_1000"))
                #plotResSeam2(res3Sp, "topright", cls = c("B","alpha","Phi","Mm"))
                #plotResSeam2(res3Sp, "topright", cls = c("EL","ER"))
                #plotResSeam2(res3Sp, "topright", cls = c("B10","EL_1000","ER_1000"))
                res3Sp$time <- res3Sp$time 
                xE3p <- tail(res3Sp,1)
                
                res3S$decRc <- res3Sc$decR
                res3S$decRp <- res3Sp$decR
                res3S$Mmc <- res3Sc$Mm
                res3S$Mmp <- res3Sp$Mm
                res3S$alphac <- res3Sc$alpha
                res3S$alphap <- res3Sp$alpha
                res3S$respRc <- with(res3Sc, resp * decR/(decR+decL))  
                res3S$respRp <- with(res3Sp, resp * decR/(decR+decL))
                res3S$scen <- scen
                res3S
            })
    resScen <- do.call( rbind, resAll)
    resScen$timeDay <- resScen$time*365
    
    endTimeDay <- 200
    dsp <- melt( subset(resScen, timeDay < endTimeDay), id=c("timeDay","scen"), measure.vars=c("decRc","decRp"),variable_name="Treatment")
    levels(dsp$Treatment) <- c("No Litter input","Litter input pulse")
    #dsp <- melt( subset(resScen, time < endTime), id=c("timeDay","scen"), measure.vars=c("L"),variable_name="Treatment")
    #dsp <- melt( subset(resScen, time < endTime), id=c("timeDay","scen"), measure.vars=c("B"),variable_name="Treatment")
    #dsp <- melt( subset(resScen, time < endTime), id=c("timeDay","scen"), measure.vars=c("limE1","limE2"),variable_name="Treatment")
    #dsp <- melt( subset(resScen, time < endTime), id=c("timeDay","scen"), measure.vars=c("ER","EL"),variable_name="Treatment")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    p3p <- ggplot( dsp, aes(x=timeDay, y=value, lty=Treatment, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (day)")+ ylab("R Depymerization (gC/m2/yr)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3p + colScale + theme( plot.margin = unit( c(0,0,0,0) , "in" ) )
    
    if (isPaperBGC){
        twWin(width=3.3, height=2.0, pointsize=9, pdf="soilPaper14/fig/PrimingDec.pdf")
        print( p3p + colScale + #theme(legend.position = c(0.9,1.05), legend.justification=c(1,1)) +
        theme( plot.margin = unit( c(0,0,0,0) , "in" ) ) )
        dev.off()
    }
    
    
    dsp <- melt( subset(resScen, timeDay < endTimeDay), id=c("timeDay","scen"), measure.vars=c("Mmc","Mmp"),variable_name="Treatment")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    levels(dsp$Treatment) <- c("No Litter input","Litter input pulse")
    dsp$value[dsp$value < 0.5] <- NA
    p3m <- ggplot( dsp, aes(x=timeDay, y=value, lty=Treatment, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (day)")+ ylab("N -Mineralization (gN/m2/yr)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3m + colScale

    if (isPaperBGC){
        twWin(width=3.3, height=3.0, pointsize=9, pdf="soilPaper14/fig/PrimingMin.pdf")
        print(p3m + colScale + theme(legend.position = c(0.9,0.53), legend.justification=c(1,1)) + theme( plot.margin = unit( c(0,0,0,0) , "in" ) ) +
        theme( plot.margin = unit( c(0,0,0,0) , "in" ) ) )
        dev.off()
    }
    
    dsp <- melt( subset(resScen, timeDay < endTimeDay), id=c("timeDay","scen"), measure.vars=c("alphac","alphap"),variable_name="Treatment")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    levels(dsp$Treatment) <- c("No Litter input","Litter input pulse")
    p3a <- ggplot( dsp, aes(x=timeDay, y=value, lty=Treatment, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (day)")+ ylab("Allocation to R: alpha (gN/m2/yr)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3a + colScale
    
    
}

simBareSoil <- function(
### Simulated decrease of input to 1/100 during years 10-210
){
    #scen <- "Revenue"
    resAll <- lapply( c("Revenue","Fixed"), function(scen){
                parmsInit <- parmsScen[[scen]]
                # spinup run
                times <- seq(0,500, length.out=101)
                #times <- seq(0,10000, length.out=101)
                res <- res1 <- as.data.frame(lsoda( x0, times, derivSeam2, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivSeam2, parms=within(parms0, useFixedAlloc<-TRUE) ))
                xE <- unlist(tail(res,1))
                plotResSeam1(res, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
                
                # 10 yr steady state
                t1S <- 10
                res <- res1S <- as.data.frame(lsoda( xE[1:length(x0)+1], 1:t1S, derivSeam2, parms=parmsInit))
                
                # 30 yr decreased C input    
                t2I = 200
                fInputInc = 1.2
                parmsC2 <- within(parmsInit, { # double as much C in 
                            iL <- iL/100
                            #cnIL <- cnIL*fInputInc
                            #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
                        }) 
                times <- seq(0,t2I, length.out=101)
                res <- res2I <- as.data.frame(lsoda( xE[1:length(x0)+1], times, derivSeam2, parms=parmsC2))
                res2I$time <- res2I$time +t1S 
                xE2 <- unlist(tail(res2I,1))
                plotResSeam1(res, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
                #trace(derivSeam2, recover)        #untrace(derivSeam2)
                tmp <- derivSeam2(0, xE2[1:length(x0)+1], parmsC2)
                
                # 30 yr normal input
                t3S <- 50
                times <- seq(0,t3S, length.out=101)
                res <- res3S <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivSeam2, parms=parmsInit))
                res3S$time <- res3S$time +t1S + t2I
                xE3 <- tail(res3S,1)
                plotResSeam1(res, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
                
                resc <- rbind( res1S, res2I[-1,], res3S[-1,])
                resc$scen <- scen
                plotResSeam1(resc, "topright", cls = c("B10","respO","Mm","Rr","Lr","alpha100"))
                resc$tvrR <- 1/(parmsInit$kR * resc$limER)
                resc
            })
    resScen <- do.call( rbind, resAll)
    resScen$Allocation <- factor(resScen$scen, levels=c("Fixed","Match","Revenue"))
    
    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("R","L"),variable_name="Pool")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    #names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area() 
    
    p3 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3 + colScale
    
    
    p3b <- ggplot( resScen, aes(x=time, y=tvrR, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Turnover time R (yr)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            #ylim(c(0,800)) +
            theme()                
    p3b + colScale
    
}




#isPaperBGC <- TRUE
if (isPaperBGC){
    simfCNGraph()
    simInitSteady()
    simCO2Increase()
    simPriming()
}

