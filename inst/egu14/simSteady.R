# simulating CO2 increase and bare soil, based on modEezy5
# based on 
library(ggplot2)
library(reshape)  # melt
library(RColorBrewer) # brewer.pal

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
        ,useFixedAlloc=FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
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
        flexible = parms0
        ,fixed = within(parms0, {isAlphaFix <- TRUE})
        ,match = within(parms0, {isAlphaMatch <- TRUE})
)

# use the same colors in every graph - color-blind friendly
myColors <- brewer.pal(5,"Dark2")
names(myColors) <- names(parmsScen)
colScale <- scale_colour_manual(name = "Allocation",values = myColors)

simInitSteady <- function(){
    scen <- "flexible"
    resAll <- lapply( c("flexible","fixed","match"), function(scen){
                parmsInit <- parmsScen[[scen]]
                times <- seq(0,100, length.out=101)
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
            theme_bw(base_size=10) +
            theme()                
    p1+ colScale
    
}
                
simCO2Increase <- function(){
     #scen <- "flexible"
     resAll <- lapply( c("flexible","fixed"), function(scen){
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
            
            # 30 yr double C input    
            t2I = 50
            fInputInc = 1.2
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
            t3S <- 50
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
    dsp$Allocation <- factor(dsp$scen, levels=c("fixed","match","flexible"))
    #names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area() 
    
    p2 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") +
            theme_bw(base_size=10) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p2 + colScale
}
    
simBareSoil <- function(){
    #scen <- "flexible"
    resAll <- lapply( c("flexible","fixed"), function(scen){
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
                resc
            })
    resScen <- do.call( rbind, resAll)
    
    dsp <- melt(resScen, id=c("time","scen"), measure.vars=c("S1","S2"),variable_name="Pool")
    dsp$Allocation <- factor(dsp$scen, levels=c("fixed","match","flexible"))
    #names(dsp)[names(dsp)=="scen"] <- "Allocation"
    #p1 <- ggplot( dsp, aes(x=time, y=value, fill=Pool, lty=scen)) + geom_area() 
    
    p3 <- ggplot( dsp, aes(x=time, y=value, lty=Pool, col=Allocation)) + geom_line(size=1) + 
            xlab("Time (yr)")+ ylab("Carbon stock (gC/m2)") +
            theme_bw(base_size=10) +
            #scale_colour_discrete(drop=TRUE,limits = levels(dsp$Allocation)) +
            theme()                
    p3 + colScale
}





