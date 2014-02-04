#library(deSolve)
# extension to have really two enzyme pools, only synthesis is split
# here with S2 pool having also a CN ratio -> different calculation of alpha

# problems with switching alpha between two values

# mg C and N, days
x0 <- x0Orig <- c(
        B = 150       ##<< microbial biomass 
        #,DOM = 0    ##<< dissolved organic matter
        ,E1  = 0.1  ##<< total enzyme pool
        ,E2  = 0.1  ##<< total enzyme pool
        #,S1 = 500    ##<< N rich substrate
        #,S2 = 500    ##<< N poor substrate
        ,I =  0      ##<< inorganic pool
)
x <- x0



parms0 <- list(
        cnB = 7.16
        ,cnE = 7.16     #tom 3
        ,cnS1 = 5
        ,cnS2 = 100
        ,kN = 0.05   ##<< (per day) enzyme turnover
        #,kS0 = 1     ##<< substrate decomposition rate (kS will be decomposition flux)
        ,kS1 = 5e-3  ##<< substrate decomposition rate N-rich (here simulate large N stock)
        #,kS2 = 10e-3 ##<< substrate decomposition rate N-poor (decomposition faster, advantage for C)
        ,kS2 = 5e-2 ##<< substrate decomposition rate N-poor (decomposition faster, advantage for C)
        ,aE = 0.05   ##<< C-uptake allocated to enzymes
        ,K = 0.3    ##<< enzyme half-saturation constant
        ,m = 0.01    ##<< maintenance respiration rate
        ,tau = 0.012    ##<< biomass turnover rate
        ,eps = 0.5      ##<< carbon use efficiency
        ,useAlpha0 = FALSE  ##<< flag to indicate equation without maintenance
)
parms0 <- within(parms0,{
 #kS1 <- kS2 <- kS <- kS0*500
            kS1S <- kS1*500
            kS2S <- kS2*100
            K1 <- K2 <- K
 eps1 <- eps2 <- eps
        })

parms <- parms0
x <- x0


.calcAlphaM <- function(E, respMaint, parms){
    # calculating alpha with maintenance
    asqrt <- E^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 2*E^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + E^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 2 *E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 4*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms $cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 *parms$kS*respMaint + E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2  + 2*E*parms$K*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 4*E*parms$K*parms$cnB^2 *parms$cnS1*parms$cnS2*parms$kS^2 + 2*E*parms$K*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 +  6*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 8 *E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 6*E*parms$ K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 8*E*parms$K* parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^ 2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 - 2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + parms$K^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1^2 *parms$cnS2*parms$eps*parms$kS*respMaint + 4*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2
    alpha <- alpha_raw <- if( asqrt < 0 ){
                #stop("encountered square root < 0")
                alpha <- 1  # take all from N-rich
            }else {
                alpha <- -(-E^2*parms$cnB*parms$cnS1*parms$kS - E^2*parms$cnB*parms$cnS2*parms$kS + 2*E ^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS - E^2*parms$cnS1*parms$cnS2*parms$ eps*respMaint + E*parms$K*parms$cnB*parms$cnS1*parms$kS - E*parms$K*parms$cnB*parms$cnS2*parms$kS + 
                            sqrt( asqrt)*abs(E))/(2*E^2*parms$cnB*parms$cnS1*parms$kS + 2*E^2* parms$cnB*parms$cnS2*parms$kS 
                            - 4*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS + 2*E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint)
            }
    #alpha <- alpha_raw2 <- (E^2*parms$cnB*parms$cnS1*parms$kS + E^2*parms$cnB*parms$cnS2*parms$kS - 2*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS + E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint - E*parms$K*parms$cnB*parms$cnS1*parms$kS + E*parms$K*parms$cnB*parms$cnS2*parms$kS + sqrt(E^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 2*E^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + E^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 4*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 4*E*parms$K*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + 2*E*parms$K*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 8*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 - 2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + parms$K^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint + 4*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$cnS1*parms$kS + 2*E^2*parms$cnB*parms$cnS2*parms$kS - 4*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS +2*E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint)
    if( alpha < 0 ) alpha = 0
    if( alpha > 1 ) alpha = 1
    alpha
}

.calcAlphaNoM <- function(E, parms){
    # calculating alpha without maintenance
    alpha0 <- (E*parms$cnB*parms$cnS1 + E*parms$cnB*parms$cnS2 - 2*E*parms$cnS1*parms$cnS2*parms$eps - parms$K*parms$cnB*parms$cnS1 + parms$K*parms$cnB*parms$cnS2 - sqrt(E^2*parms$cnB^2*parms$cnS1^2 + 2*E^2*parms$cnB^2*parms$cnS1*parms$cnS2 + E^2*parms$cnB^2*parms$cnS2^2 - 4*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps- 4*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 2*E*parms$K*parms$cnB^2*parms$cnS1^2 + 4*E*parms$K*parms$cnB^2*parms$cnS1*parms$cnS2 + 2*E*parms$K*parms$cnB^2*parms$cnS2^2 - 8*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps - 8*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 8*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + parms$K^2*parms$cnB^2*parms$cnS1^2 - 2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2 + parms$K^2*parms$cnB^2*parms$cnS2^2))/(2*E*parms$cnB*parms$cnS1 + 2*E*parms$cnB*parms$cnS2 - 4*E*parms$cnS1*parms$cnS2*parms$eps)
    if( alpha < 0 ) alpha = 0
    if( alpha > 1 ) alpha = 1
    alpha
}


.tmp.fa <- function(E){
    fAlpha <- function(E){
        (E^2*parms$cnB*parms$cnS1*parms$kS + E^2*parms$cnB*parms$cnS2*parms$kS - 2*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS + E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint - E*parms$K*parms$cnB*parms$cnS1*parms$kS + E*parms$K*parms$cnB*parms$cnS2*parms$kS + sqrt(E^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 2*E^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + E^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 4*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 4*E*parms$K*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + 2*E*parms$K*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 8*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 - 2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + parms$K^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint + 4*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$cnS1*parms$kS + 2*E^2*parms$cnB*parms$cnS2*parms$kS - 4*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS +2*E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint)
        #-(-E^2*parms$cnB*parms$cnS1*parms$kS - E^2*parms$cnB*parms$cnS2*parms$kS + 2*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS - E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint + E*parms$K*parms$cnB*parms$cnS1*parms$kS - E*parms$K*parms$cnB*parms$cnS2*parms$kS + sqrt(E^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 2*E^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + E^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 4*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 4*E*parms$K*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + 2*E*parms$K*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 8*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 - 2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + parms$K^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint + 4*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$cnS1*parms$kS + 2*E^2*parms$cnB*parms$cnS2*parms$kS - 4*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS+ 2*E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint)
    }
    Es <- 1:500
    Alphas <- fAlpha(Es)
    plot( Alphas ~ Es)
}

derivEezy2b <- function(t,x,parms){
    E <- x["E1"]+x["E2"]
    # x["E1"]/E
    respMaint <- parms$m * x["B"]
    #respMaint <- 0
    # including maintenance
    #alpha <- -(-E^2*parms$cnB*parms$kS - 2*E^2*parms$cnS1*parms$eps*parms$kS - E^2*parms$cnS1*parms$eps*respMaint - E*parms$K*parms$cnB*parms$kS + sqrt(E^2*parms$cnB^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$kS - 4*E^2*parms$cnS1*parms$eps*parms$kS + 2*E^2*parms$cnS1*parms$eps*respMaint)    
    #x[c("B","E1","E2","S1","S2")] <- pmax(0,x[c("B","E1","E2","S1","S2")])
    #alpha2 <- 1/2 + sqrt(-E^2*(E + 2*parms$K)*(2*parms$eps*parms$kS - respMaint)*(-2*E*parms$eps*parms$kS + E*respMaint + 2*parms$K*respMaint))/(2*E^2*(2*parms$eps*parms$kS - respMaint))
    alpha <- calcAlphaMS(E,respMaint, parms$kS1S, parms$kS2S, parms)
    E1 <- x["E1"]
    E2 <- x["E2"]
    decC1 <- parms$kS1S *  E1 / (parms$K + E1)
    decC2 <- parms$kS2S *  E2 / (parms$K + E2)
    parms$eps * (decC1+decC2)/(decC1/parms$cnS1 + decC2/parms$cnS2)
    parms$eps * (decC1+decC2-respMaint)/(decC1/parms$cnS1 + decC2/parms$cnS2)
    uC <- decC1 + decC2
    uN <- decC1/parms$cnS1 + decC2/parms$cnS2
    #(uC)* pr$eps/ (uN)         # check to be cnB = 7.16
    minUCE <- min(uC, uN*parms$cnE)  # C for enzyme biosynthesis, if necessary take microbial material for associated resp
    isLimEN <- if(uN*parms$cnE < uC) 1 else 0
    synE <- parms$aE * minUCE
    respSynE <- (1-parms$eps)/parms$eps * synE
    ucmC <- uC-synE-respSynE-respMaint
    ucmN <- (uN-synE/parms$cnE)*parms$cnB/parms$eps
    isLimN <- if(ucmN < ucmC) 1 else 0
    minUCM <- min( ucmC, ucmN ) # C flux available for biomass synthesis
    if( minUCM > 0){    
        synB <- parms$eps * minUCM  
        respSynB <- (1-parms$eps)/parms$eps * synB
        edecay <- 0
    }else{
        #print("negative C balance; recover"); recover()
        synB <- minUCM      # negative: endogeneous decay, releases N
        respSynB <- 0
    }
    tvrB <- parms$tau*x["B"]
    respO <- uC - (synE+synB+respSynE+respSynB+respMaint)
    Mm <- uN - (synE/parms$cnE + synB/parms$cnB)
         
    #if( synB < tvrB | parms$useAlpha0){
    if( FALSE ){
        # provide CN/ratio for enzyme production before maintenance
        alpha0 <- (E*parms$cnB*parms$cnS1 + E*parms$cnB*parms$cnS2 - 2*E*parms$cnS1*parms$cnS2*parms$eps - parms$K*parms$cnB*parms$cnS1 + parms$K*parms$cnB*parms$cnS2 - sqrt(E^2*parms$cnB^2*parms$cnS1^2 + 2*E^2*parms$cnB^2*parms$cnS1*parms$cnS2 + E^2*parms$cnB^2*parms$cnS2^2 - 4*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps- 4*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 2*E*parms$K*parms$cnB^2*parms$cnS1^2 + 4*E*parms$K*parms$cnB^2*parms$cnS1*parms$cnS2 + 2*E*parms$K*parms$cnB^2*parms$cnS2^2 - 8*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps - 8*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 8*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + parms$K^2*parms$cnB^2*parms$cnS1^2 - 2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2 + parms$K^2*parms$cnB^2*parms$cnS2^2))/(2*E*parms$cnB*parms$cnS1 + 2*E*parms$cnB*parms$cnS2 - 4*E*parms$cnS1*parms$cnS2*parms$eps)
        alpha <- alpha0        
    } 
    dB <- synB - tvrB
    dE1 <- + alpha*synE  - parms$kN*x["E1"]
    dE2 <- + (1-alpha)*synE  - parms$kN*x["E2"]
    #dS1 <- -decC1
    #dS2 <- -decC2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    
    if( abs(diff( cBal <- c(uC=uC, usage=resp+synB+synE ))) > 1e-6 ) stop("C-balance error")  # check the same
    if( abs(diff( cBal <- c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + Mm ) )) > 1e-6) stop("N-balance error") # check the same
    
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dI)),names=c("dB","dE1","dE2","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    #, alpha2 = as.numeric(alpha2)
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm)), alpha=as.numeric(alpha), isLimEN=isLimEN, isLimN=isLimN, decC1=as.numeric(decC1), decC2=as.numeric(decC2) )
}


.tmp.f <- function(){
    times <- seq(0,500, length.out=101)
    derivEezy2b(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    parmsMsmall <- within(parms0, m <- 1e-3)   # small maintenance
    #res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy2b, parms=parmsM0))
    res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy2b, parms=parms0))
    #res <- res1 <- as.data.frame(lsoda( x0, times, derivEezy2b, parms=parms0,atol = 1e-4))
    
    .tmp.f <- function(){
        #trace(derivEezy2b, recover)  # untrace(derivEezy2b)
        derivEezy2b(0, tail(res,1)[2:5], parms0)
    }
    
    res$B100 <- res$B/100
    res$E1_10 <- res$E1/10
    res$E2_10 <- res$E2/10
    res$respO_10 <- res$respO/10
    cls <- c("B100","E1","E2","respO","Mm")
    #cls <- c("B1000","E1_10","E2_10","respO_10","Mm")
    #cls <- c("B1000","E1_10","E2_10","respO","Mm","alpha")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:10, col=1:10 )  #, ylim=c(0,30))
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
}

.tmp.fCNGraph <- function(){
    fdEezy <- derivEezy2b
    #fdEezy <- derivEezy4bc
    
    
    times <- seq(0,500, length.out=101)
    cnS1s <- seq(5,25,by=.2)
    
    .tmp.f <- function(){
        # no decomposition of S2        
        parmsNoS2 <- within(parms0, kS2 <- 0)
        cnAll <- cnS1s
        resl <- lapply(cnS1s, function(cnS1i){
                    parms <- within(parmsNoS2, cnS1 <- cnS1i)
                    res <- as.data.frame(lsoda( x0, times, fdEezy, parms=parms))
                    tail(res,1)
                })
    }
    
    #parms0 <- within(parms0, useFixedAlloc <- TRUE )
    #cnS1i <- cnS1s[1]
    #cnS1i <- cnS1s[length(cnS1s)]
    cnAll <- 1000/(500/cnS1s)
    resl <- lapply(cnS1s, function(cnS1i){
                parms <- within(parms0, cnS1 <- cnS1i)
                res <- as.data.frame(lsoda( x0, times, fdEezy, parms=parms))
                tail(res,1)
            })
    .tmp.f <- function(){
        #trace(derivEezy, recover)   # untrace(derivEezy)
        derivEezy( 0, tail(res,1)[2:6], parms)
    }
    
    resE1 <- do.call( rbind, resl)
    resE <- cbind( cnS1 = cnS1s, cnAll=cnAll, resE1 )
    resE$B1000 <- resE$B/1000
    resE$E1_10 <- resE$E1/10
    resE$alpha_10 <- resE$alpha/10
    cls <- c("B1000","respO","Mm","alpha_10")
    #cls <- c("B1000","E10","respO","Mm")
    matplot( cnS1s, resE[,cls], type="l", xlab="CN S1 (CN S2=100)", ylab="", ylim=c(0,0.3))
    legend("topleft", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)

    #install.packages("dplyr")
    iMax <- which.max(resE$B)
    resE[iMax+(-2:2),]
}






