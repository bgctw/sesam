# simulating CO2 increase and bare soil, based on modEezy5
# based on 
library(ggplot2)
library(reshape)  # melt
library(RColorBrewer) # brewer.pal

baseFontSize <- 10  # pubs
baseFontSize <- 16  # presentations

# gC/m2 and gN/m2, /yr
parms0 <- list(
        cnB = 7.16
        ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
        ,cnIS1 = 4.5      ##<< between micr and enzyme signal
        ,cnIS2 = 30 ##<< N poor substrate
        ,kN= 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
        ,kNB = 0.8     ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of S1)
        #,kS1 = 1/(50)        ##<< 1/(x years) 
        ,kS1 = 1/(10)        ##<< 1/(x years)       # to demonstrate changes on short time scale
        ,kS2 = 1/(1)        ##<< 1/(x years) 
        ,aE = 0.001*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
        ,km = 0.3     ##<< enzyme half-saturation constant
        #,km = 0.03     ##<< enzyme half-saturation constant
        ,m = 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
        ,tau = 1/60*365  ##<< biomass turnover rate (12 days)
        ,eps = 0.5      ##<< carbon use efficiency
        ,epsTvr = 0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
        ,iS1 = 0        ##<< input modelled explicitely
        ,iS2 = 300         # g/m2 input per year (half NPP)
        #,plantNUp = 300/70*1/4  # plant N uptake balancing N inputs
        ,plantNUp = 0
        ,useFixedAlloc=FALSE    ##<< set to true to use Fixed enzyme allocation (alpha = 0.5)
)
parms0 <- within(parms0,{
            km1 <- km2 <- km
            eps1 <- eps2 <- eps
            cnE1 <- cnE2 <- cnE 
            kN1 <- kN2 <- kN
        })
parms <- parms0

x0 <- x0Orig <- c( #aE = 0.001*365
        B = 20                     ##<< microbial biomass 
        ,E1  = 1.5*parms0$km                  ##<< total enzyme pool
        ,E2  = 1.5*parms0$km                   ##<< total enzyme pool
        ,S1 = 490                 ##<< N rich substrate
        ,SN1 = 490/parms0$cnIS1    ##<< N rich substrate N pool
        ,S2 = 380                  ##<< N poor substrate
        ,SN2 = 380/parms0$cnIS2    ##<< N poor substrate N pool
        ,I =  0                    ##<< inorganic pool
)
x <- x0


parmsScen <- list(
        Revenue = parms0
        ,Fixed = within(parms0, {isAlphaFix <- TRUE})
        ,Match = within(parms0, {isAlphaMatch <- TRUE})
)

# use the same colors in every graph - color-blind friendly
myColors <- brewer.pal(5,"Dark2")
.simpleCap <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2),  sep = "")
}
#names(myColors) <- .simpleCap(names(parmsScen))
names(myColors) <- names(parmsScen)
colScale <- scale_colour_manual(name = "Allocation",values = myColors)


simfCNGraph <- function(
    ### fixing Substrate availability and varying N amount in S2
){
    parms0F <- within(parms0,{    isFixedS <- TRUE; kNB = 0     })
    parmsScenF <- list(
            Revenue = parms0F
            ,Fixed = within(parms0F, {isAlphaFix <- TRUE})
            ,Match = within(parms0F, {isAlphaMatch <- TRUE})
    )
    
    cnS1 <- 6.8
    cnS2 <- 30
    x0N <- c( #aE = 0.001*365
            B = 20                     ##<< microbial biomass 
            ,E1  = 1.5*parms0$km                  ##<< total enzyme pool
            ,E2  = 1.5*parms0$km                   ##<< total enzyme pool
            ,S1 = 800                 ##<< N rich substrate
            ,SN1 = 800/cnS1    ##<< N rich substrate N pool
            ,S2 = 400                  ##<< N poor substrate
            ,SN2 = 400/cnS2    ##<< N poor substrate N pool
            ,I =  0                    ##<< inorganic pool
    )
    
    cnS2s <- seq( 18,42,by=1)
    cnS2 <- 23
    resLC <- lapply( cnS2s, function(cnS2){ 
                x0N["SN2"] <- x0N["S2"]/cnS2
                times <- seq(0,10, length.out=101)
                #times <- seq(0,10000, length.out=101)
                #scen <- "Match"
                resL <- lapply( names(parmsScenF), function(scen){
                #resL <- lapply( "Revenue", function(scen){
                            parmsInit <- parmsScenF[[scen]]
                            res <- res1 <- as.data.frame(lsoda( x0N, times, derivEezy5, parms=parmsInit))
                            xE <- unlist(tail(res,1))
                            #plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                            #trace(derivEezy5, recover) #untrace(derivEezy5)
                            #tmp <- derivEezy5(0, xE[1:length(x0)+1], parmsInit)
                            xE
                        })
                resScen <- cbind( data.frame(scen=names(parmsScenF), do.call( rbind, resL )))
            })
    resC <- resCAll <- do.call( rbind, resLC )
    resC <- subset(resCAll, cnS2 <= 50)
    resC$MmImbMon <- resC$MmImb/12
    resC$respOMon <- resC$respO/12

    #p8a <- ggplot(subset(resC, scen=="Revenue"), aes(x=cnS2, y=isLimN)) + geom_line()
    #p8a
    resCRev <- subset(resC, scen=="Revenue")
    cnTER <- resCRev$cnS2[ which.min( abs(1-resCRev$isLimN) ) ]
    
    dsp <- melt(resC, id=c("cnS2","scen"), measure.vars=c("MmImbMon","respOMon","alpha","B"), variable_name="Measure")
    #levels(dsp$Measure) <- c("Allocation_Ratio (alpha)","Biomass (gC/m^2)","Mineralization[Imb] (gN/m^2/yr)","OverflowRespiration (gC/m2/yr)")
    levels(dsp$Measure) <- c("Allocation to E_SOM (alpha)","Biomass (gC/m2)","Mineralization Imb (gN/m2/month)","OverflowRespiration (gC/m2/month)")
    p8 <- ggplot( dsp, aes(x=cnS2, y=value, col=scen) ) + geom_line(size=1) + 
            #facet_grid(Measure~., scales = "free", labeller= label_parsed ) + 
            facet_wrap(~Measure, scales = "free_y" ) +
            scale_x_continuous('C/N ratio of Lit') + 
            #geom_vline(aes(xintercept=cnTER), colour="#990000", linetype="dashed") +
            theme_bw(base_size=baseFontSize) +
            theme(axis.title.y = element_blank()) +
            theme()
    #twWin(6)
    p8 + colScale
    

    dsp <- melt(resC, id=c("cnS2","scen"), measure.vars=c("alpha","alphaC","alphaN","rNLim", "pCLim", "pNLim"), variable_name="Measure")
    #levels(dsp$Measure) <- c("Allocation_Ratio (alpha)","Biomass (gC/m^2)","Mineralization[Imb] (gN/m^2/yr)","OverflowRespiration (gC/m2/yr)")
    #levels(dsp$Measure) <- c("Allocation to E_SOM (alpha)","Biomass (gC/m2)","Mineralization Imb (gN/m2/month)","OverflowRespiration (gC/m2/month)")
    p8a <- ggplot( dsp, aes(x=cnS2, y=value, col=scen) ) + geom_line(size=1) + 
            #facet_grid(Measure~., scales = "free", labeller= label_parsed ) + 
            facet_wrap(~Measure, scales = "free_y" ) +
            scale_x_continuous('C/N ratio of Lit') + 
            geom_vline(aes(xintercept=cnTER), colour="#990000", linetype="dashed") +
            theme_bw(base_size=baseFontSize) +
            theme(axis.title.y = element_blank()) +
            theme()
    #twWin(6)
    p8a 
    
    
    dspC <- melt(subset(resC, scen=="Revenue"), id=c("cnS2"), measure.vars=c("effE1C","effE2C"), variable_name="Pool"); dspC$Element = "C"
    levels(dspC$Pool) <- c("SOM","Lit")
    dspN <- melt(subset(resC, scen=="Revenue"), id=c("cnS2"), measure.vars=c("effE1N","effE2N"), variable_name="Pool"); dspN$Element = "N"
    levels(dspN$Pool) <- c("SOM","Lit")
    dsp <- rbind(dspC, dspN)
    #levels(dsp$Measure) <- c("Allocation_Ratio (alpha)","Biomass (gC/m^2)","Mineralization[Imb] (gN/m^2/yr)","OverflowRespiration (gC/m2/yr)")
    #levels(dsp$Measure) <- c("Allocation ratio (alpha)","Biomass (gC/m2)","Mineralization Imb (gN/m2/month)","OverflowRespiration (gC/m2/month)")
    p8b <- ggplot( dsp, aes(x=cnS2, y=value, col=Pool) ) + geom_line(size=1) + 
            #facet_grid(Measure~., scales = "free", labeller= label_parsed ) + 
            facet_wrap(~Element, scales = "free_y" ) +
            scale_x_continuous('C/N ratio of Lit') +
            geom_vline(aes(xintercept=cnTER), colour="#990000", linetype="dashed") +
            theme_bw(base_size=baseFontSize) +
            #theme(axis.title.y = element_blank()) +
            ylab("Revenue")
            theme()
    #twWin(6)
    p8b 
    
    
    
}


simInitSteady <- function(
    ### inspect approaching a steady state (or breakdown of biomass)
){
    scen <- "Revenue"
    resAll <- lapply( c("Revenue","Fixed","Match"), function(scen){
                parmsInit <- parmsScen[[scen]]
                times <- seq(0,70, length.out=101)
                #times <- seq(0,10000, length.out=101)
                res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivEezy5, parms=within(parms0, useFixedAlloc<-TRUE) ))
                xE <- unlist(tail(res,1))
                plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                res$scen <- scen
                res
            })
    resScen <- do.call( rbind, resAll)
    
    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("S1","S2"),variable_name="Pool")
    names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area() 
    
    p1 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") +
            theme_bw(base_size=baseFontSize) +
            ylim(c(300,600))
            theme()                
    p1+ colScale
    
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
            # spinup run
            times <- seq(0,500, length.out=101)
            #times <- seq(0,10000, length.out=101)
            res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parmsInit))
            #res <- res1f <- as.data.frame(lsoda( x0, times, derivEezy5, parms=within(parms0, useFixedAlloc<-TRUE) ))
            xE <- unlist(tail(res,1))
            plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
            
            # 10 yr steady state
            res <- res1S <- as.data.frame(lsoda( xE[1:length(x0)+1], 1:t1S, derivEezy5, parms=parmsInit))
            
            # 30 yr double C input    
            parmsC2 <- within(parmsInit, { # double as much C in 
                        iS2 <- iS2*fInputInc
                        cnIS2 <- cnIS2*fInputInc
                        #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
                    }) 
            times <- seq(0,t2I, length.out=101)
            res <- res2I <- as.data.frame(lsoda( xE[1:length(x0)+1], times, derivEezy5, parms=parmsC2))
            res2I$time <- res2I$time +t1S 
            xE2 <- unlist(tail(res2I,1))
            plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
            #trace(derivEezy5, recover)        #untrace(derivEezy5)
            tmp <- derivEezy5(0, xE2[1:length(x0)+1], parmsC2)
            
            # 30 yr normal input
            times <- seq(0,t3S, length.out=101)
            res <- res3S <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivEezy5, parms=parmsInit))
            res3S$time <- res3S$time +t1S + t2I
            xE3 <- tail(res3S,1)
            plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
            
            resc <- rbind( res1S, res2I[-1,], res3S[-1,])
            resc$scen <- scen
            plotRes(resc, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
            resc
    })
    resScen <- do.call( rbind, resAll)

    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("S1","S2"),variable_name="Pool")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    #names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area()

    p2 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p2 + colScale
    
    dsCN <- data.frame( time=1:(t1S+t2I+t3S), iS2=parmsScen[["Revenue"]]$iS2 )
    dsCN[ dsCN$time %in% (t1S+1):(t1S+t2I),"iS2"] <- parmsScen[["Revenue"]]$iS2 * fInputInc

    p2b <- ggplot( dsCN, aes(x=time, y=iS2)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Input") + ylim(c(0,400)) +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()
    p2b
    
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
    #scen <- "Revenue"
    resAll <- lapply( c("Revenue","Fixed"), function(scen){
                parmsInit <- parmsScen[[scen]]
                
                # 500 yr decreased C input    
                t2I = 500
                fInputInc = 1.2
                parmsC2 <- within(parmsInit, { # double as much C in 
                            iS2 <- iS2/100
                            #cnIS2 <- cnIS2*fInputInc
                            #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
                        }) 
                times <- c(seq(0,t2I, length.out=101), 1:10+t2I)
                res <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parmsC2))
                res2I <- tail(res,10)   # last 10 years
                xE2 <- unlist(tail(res2I,1))
                #plotRes(res2I, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                #plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                #trace(derivEezy5, recover)        #untrace(derivEezy5)
                #tmp <- derivEezy5(0, xE2[1:length(x0)+1], parmsC2)
                
                # 5 yr control of continued decreased input 
                #t3S <- 5
                t3S <- 50
                times <- seq(0,t3S, length.out=201)
                res <- res3Sc <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivEezy5, parms=parmsC2))
                res3Sc$time <- res3Sc$time 
                xE3c <- tail(res3Sc,1)
                
                # 5 yr of continued decreased input but with initial pulse of S2
                x0P <- xE2; x0P["S2"] <- x0P["S2"] + 50 
                res <- res3S <- res3Sp <- as.data.frame(lsoda( x0P[1:length(x0)+1], times, derivEezy5, parms=parmsC2))
                res3Sp$time <- res3Sp$time 
                xE3p <- tail(res3Sp,1)
                
                res3S$decC1c <- res3Sc$decC1
                res3S$decC1p <- res3Sp$decC1
                res3S$Mmc <- res3Sc$Mm
                res3S$Mmp <- res3Sp$Mm
                res3S$respS1c <- with(res3Sc, resp * decC1/(decC1+decC2))  
                res3S$respS1p <- with(res3Sp, resp * decC1/(decC1+decC2))
                res3S$scen <- scen
                res3S
            })
    resScen <- do.call( rbind, resAll)
    
    dsp <- melt( subset(resScen, time < 8), id=c("time","scen"), measure.vars=c("decC1c","decC1p"),variable_name="Treatment")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    levels(dsp$Treatment) <- c("No Litter input","Litter input pulse")
    p3p <- ggplot( dsp, aes(x=time, y=value, lty=Treatment, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("SOM Decompos. (gC/m2/yr)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3p + colScale
    
    dsp <- melt( subset(resScen, time < 8), id=c("time","scen"), measure.vars=c("Mmc","Mmp"),variable_name="Treatment")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    levels(dsp$Treatment) <- c("No Litter input","Litter input pulse")
    p3p <- ggplot( dsp, aes(x=time, y=value, lty=Treatment, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Mineralization (gN/m2/yr)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3p + colScale
    
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
                res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy5, parms=parmsInit))
                #res <- res1f <- as.data.frame(lsoda( x0, times, derivEezy5, parms=within(parms0, useFixedAlloc<-TRUE) ))
                xE <- unlist(tail(res,1))
                plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                
                # 10 yr steady state
                t1S <- 10
                res <- res1S <- as.data.frame(lsoda( xE[1:length(x0)+1], 1:t1S, derivEezy5, parms=parmsInit))
                
                # 30 yr decreased C input    
                t2I = 200
                fInputInc = 1.2
                parmsC2 <- within(parmsInit, { # double as much C in 
                            iS2 <- iS2/100
                            #cnIS2 <- cnIS2*fInputInc
                            #plantNUp <- 300/70*4/5  # plant N uptake balancing N inputs
                        }) 
                times <- seq(0,t2I, length.out=101)
                res <- res2I <- as.data.frame(lsoda( xE[1:length(x0)+1], times, derivEezy5, parms=parmsC2))
                res2I$time <- res2I$time +t1S 
                xE2 <- unlist(tail(res2I,1))
                plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                #trace(derivEezy5, recover)        #untrace(derivEezy5)
                tmp <- derivEezy5(0, xE2[1:length(x0)+1], parmsC2)
                
                # 30 yr normal input
                t3S <- 50
                times <- seq(0,t3S, length.out=101)
                res <- res3S <- as.data.frame(lsoda( xE2[1:length(x0)+1], times, derivEezy5, parms=parmsInit))
                res3S$time <- res3S$time +t1S + t2I
                xE3 <- tail(res3S,1)
                plotRes(res, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                
                resc <- rbind( res1S, res2I[-1,], res3S[-1,])
                resc$scen <- scen
                plotRes(resc, "topright", cls = c("B10","respO","Mm","S1r","S2r","alpha100"))
                resc$tvrS1 <- 1/(parmsInit$kS1 * resc$limE1)
                resc
            })
    resScen <- do.call( rbind, resAll)
    resScen$Allocation <- factor(resScen$scen, levels=c("Fixed","Match","Revenue"))
    
    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("S1","S2"),variable_name="Pool")
    dsp$Allocation <- factor(dsp$scen, levels=c("Fixed","Match","Revenue"))
    #names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area() 
    
    p3 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3 + colScale
    
    
    p3b <- ggplot( resScen, aes(x=time, y=tvrS1, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Turnover time S1 (yr)") +
            theme_bw(base_size=baseFontSize) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            ylim(c(0,800)) +
            theme()                
    p3b + colScale
    
}




