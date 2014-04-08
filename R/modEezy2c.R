#library(deSolve)
# same as modeEezy2b but now actually tracking substrate pools


# problems with switching alpha between two values

# mg C and N, days
x0 <- x0Orig <- c(
        B = 25000       ##<< microbial biomass 
        #,DOM = 0    ##<< dissolved organic matter
        ,E1  = 0.01  ##<< total enzyme pool
        ,E2  = 0.01  ##<< total enzyme pool
        ,S1 = 500    ##<< N rich substrate
        ,S2 = 100    ##<< N poor substrate
        ,I =  0      ##<< inorganic pool
)
x <- x0



parms0 <- list(
        cnB = 7.16
        ,cnE = 7.16     #tom 3
        ,cnS1 = 5
        ,cnS2 = 100
        ,kN = 0.05   ##<< (per day) enzyme turnover
        ,kS1 = 5e-3  ##<< substrate decomposition rate N-rich (here simulate large N stock)
        #,kS2 = 10e-3 ##<< substrate decomposition rate N-poor (decomposition faster, advantage for C)
        ,kS2 = 5e-2 ##<< substrate decomposition rate N-poor (decomposition faster, advantage for C)
        #,kS0 = 0.01     ##<< substrate decomposition rate (kS will be decomposition flux)
        ,aE = 0.05   ##<< C-uptake allocated to enzymes
        ,K = 0.3    ##<< enzyme half-saturation constant
        ,m = 0.01    ##<< maintenance respiration rate
        ,tau = 0.012    ##<< biomass turnover rate
        ,eps = 0.5      ##<< carbon use efficiency
        ,useAlpha0 = FALSE  ##<< flag to indicate equation without maintenance
        ,iS1 = 1
        ,iS2 = 1
)
parms0 <- within(parms0,{
 #kS1 <- kS2 <- kS <- kS0
 K1 <- K2 <- K
 eps1 <- eps2 <- eps
        })

parms <- parms0
x <- x0


calcAlphaMS <- function(E, respMaint, decPotS1, decPotS2, parms){
    # calculating alpha with maintenance
    asqrt <- E^2*decPotS1^2*parms$cnB^2*parms$cnS2^2 - 
            2*E^2*decPotS1^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 
            E^2*decPotS1^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 2*E^2*decPotS1*decPotS2*parms$cnB^2*parms$cnS1*parms$cnS2 - 
            2*E^2*decPotS1*decPotS2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps - 2*E^2*decPotS1*decPotS2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 
            2*E^2*decPotS1*decPotS2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 2*E^2*decPotS1*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*respMaint - 
            2*E^2*decPotS1*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint + E^2*decPotS2^2*parms$cnB^2*parms$cnS1^2 - 
            2*E^2*decPotS2^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps + E^2*decPotS2^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 
            2*E^2*decPotS2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*respMaint - 2*E^2*decPotS2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint + 
            E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + 2*E*decPotS1^2*parms$K*parms$cnB^2*parms$cnS2^2 - 
            4*E*decPotS1^2*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 2*E*decPotS1^2*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 
            4*E*decPotS1*decPotS2*parms$K*parms$cnB^2*parms$cnS1*parms$cnS2 - 4*E*decPotS1*decPotS2*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps - 
            4*E*decPotS1*decPotS2*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 4*E*decPotS1*decPotS2*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 
            6*E*decPotS1*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*respMaint - 6*E*decPotS1*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint + 
            2*E*decPotS2^2*parms$K*parms$cnB^2*parms$cnS1^2 - 4*E*decPotS2^2*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps + 
            2*E*decPotS2^2*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 6*E*decPotS2*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*respMaint - 
            6*E*decPotS2*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint + 4*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + 
            decPotS1^2*parms$K^2*parms$cnB^2*parms$cnS2^2 -2*decPotS1^2*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps + 
            decPotS1^2*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 - 2*decPotS1*decPotS2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2 + 
            2*decPotS1*decPotS2*parms$K^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps + 2*decPotS1*decPotS2*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps - 
            2*decPotS1*decPotS2*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 4*decPotS1*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*respMaint - 
            4*decPotS1*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint + decPotS2^2*parms$K^2*parms$cnB^2*parms$cnS1^2- 
            2*decPotS2^2*parms$K^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps + decPotS2^2*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2 + 
            4*decPotS2*parms$K^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*respMaint - 4*decPotS2*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint + 
            4*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2
    alpha <- alpha_raw <- if( asqrt < 0 ){
                #stop("encountered square root < 0")
                alpha <- 1  # take all from N-rich
            }else {
                alpha <- -(-E^2*decPotS1*parms$cnB*parms$cnS2 + E^2*decPotS1*parms$cnS1*parms$cnS2*parms$eps - 
                            E^2*decPotS2*parms$cnB*parms$cnS1 + E^2*decPotS2*parms$cnS1*parms$cnS2*parms$eps - 
                            E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint - E*decPotS1*parms$K*parms$cnB*parms$cnS2 + 
                            E*decPotS1*parms$K*parms$cnS1*parms$cnS2*parms$eps + E*decPotS2*parms$K*parms$cnB*parms$cnS1 - 
                            E*decPotS2*parms$K*parms$cnS1*parms$cnS2*parms$eps + 
                            sqrt(asqrt)*
                            abs(E))/(2*E^2*decPotS1*parms$cnB*parms$cnS2 - 
                            2*E^2*decPotS1*parms$cnS1*parms$cnS2*parms$eps + 2*E^2*decPotS2*parms$cnB*parms$cnS1 - 2*E^2*decPotS2*parms$cnS1*parms$cnS2*parms$eps + 
                            2*E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint)
            }
    #alpha <- alpha_raw2 <- (E^2*parms$cnB*parms$cnS1*parms$kS + E^2*parms$cnB*parms$cnS2*parms$kS - 2*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS + E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint - E*parms$K*parms$cnB*parms$cnS1*parms$kS + E*parms$K*parms$cnB*parms$cnS2*parms$kS + sqrt(E^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 2*E^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + E^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 4*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$cnS1^2*parms$kS^2 + 4*E*parms$K*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + 2*E*parms$K*parms$cnB^2*parms$cnS2^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint - 8*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$cnS1^2*parms$kS^2 - 2*parms$K^2*parms$cnB^2*parms$cnS1*parms$cnS2*parms$kS^2 + parms$K^2*parms$cnB^2*parms$cnS2^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1^2*parms$cnS2*parms$eps*parms$kS*respMaint + 4*parms$K^2*parms$cnB*parms$cnS1*parms$cnS2^2*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$cnS2^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$cnS1*parms$kS + 2*E^2*parms$cnB*parms$cnS2*parms$kS - 4*E^2*parms$cnS1*parms$cnS2*parms$eps*parms$kS +2*E^2*parms$cnS1*parms$cnS2*parms$eps*respMaint)
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

derivEezy2c <- function(t,x,parms){
    E <- x["E1"]+x["E2"]
    # x["E1"]/E
    respMaint <- parms$m * x["B"]
    #respMaint <- 0
    # including maintenance
    #alpha <- -(-E^2*parms$cnB*parms$kS - 2*E^2*parms$cnS1*parms$eps*parms$kS - E^2*parms$cnS1*parms$eps*respMaint - E*parms$K*parms$cnB*parms$kS + sqrt(E^2*parms$cnB^2*parms$kS^2 - 4*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 2*E^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS^2 - 4*E^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + E^2*parms$cnS1^2*parms$eps^2*respMaint^2 + 2*E*parms$K*parms$cnB^2*parms$kS^2 - 8*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS^2 + 6*E*parms$K*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint + 8*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS^2 - 12*E*parms$K*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*E*parms$K*parms$cnS1^2*parms$eps^2*respMaint^2 + parms$K^2*parms$cnB^2*parms$kS^2 + 4*parms$K^2*parms$cnB*parms$cnS1*parms$eps*parms$kS*respMaint - 8*parms$K^2*parms$cnS1^2*parms$eps^2*parms$kS*respMaint + 4*parms$K^2*parms$cnS1^2*parms$eps^2*respMaint^2)*abs(E))/(2*E^2*parms$cnB*parms$kS - 4*E^2*parms$cnS1*parms$eps*parms$kS + 2*E^2*parms$cnS1*parms$eps*respMaint)    
    #x[c("B","E1","E2","S1","S2")] <- pmax(0,x[c("B","E1","E2","S1","S2")])
    #alpha2 <- 1/2 + sqrt(-E^2*(E + 2*parms$K)*(2*parms$eps*parms$kS - respMaint)*(-2*E*parms$eps*parms$kS + E*respMaint + 2*parms$K*respMaint))/(2*E^2*(2*parms$eps*parms$kS - respMaint))
    decPotS1 <- parms$kS1 * x["S1"]
    decPotS2 <- parms$kS2 * x["S2"]
    alpha <- calcAlphaMS(E,respMaint, decPotS1, decPotS2, parms)
    E1 <- x["E1"]
    E2 <- x["E2"]
    decC1 <- decPotS1 *  E1 / (parms$K + E1)
    decC2 <- decPotS2 *  E2 / (parms$K + E2)
    #parms$eps * (decC1+decC2)/(decC1/parms$cnS1 + decC2/parms$cnS2)
    #parms$eps * (decC1+decC2-respMaint)/(decC1/parms$cnS1 + decC2/parms$cnS2)
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
         
    dB <- synB - tvrB
    dE1 <- + alpha*synE  - parms$kN*x["E1"]
    dE2 <- + (1-alpha)*synE  - parms$kN*x["E2"]
    dS1 <- -decC1 +parms$iS1
    dS2 <- -decC2 +parms$iS2
    dI <- +Mm
    resp <- respSynE + respSynB + respMaint + respO
    
    if( abs(diff( cBal <- c(uC=uC, usage=resp+synB+synE ))) > 1e-6 ) stop("C-balance error")  # check the same
    if( abs(diff( cBal <- c(uN=uN, usage=synE/parms$cnE + synB/parms$cnB + Mm ) )) > 1e-6) stop("N-balance error") # check the same
    
    resDeriv <-structure(as.numeric(c( dB, dE1, dE2, dS1, dS2, dI)),names=c("dB","dE1","dE2","dS1","dS2","dI"))
    if( any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
    #, alpha2 = as.numeric(alpha2)
    list( resDeriv, c(respO=as.numeric(respO), Mm=as.numeric(Mm)), alpha=as.numeric(alpha), isLimEN=isLimEN, isLimN=isLimN, decC1=as.numeric(decC1), decC2=as.numeric(decC2) )
}


.tmp.f <- function(){
    times <- seq(0,500, length.out=101)
    derivEezy2c(0, x0, parms0)
    parmsM0 <- within(parms0, m <- 0)   # no maintenance
    parmsMsmall <- within(parms0, m <- 1e-3)   # small maintenance
    res <- res0 <- as.data.frame(lsoda( x0, times, derivEezy2c, parms=parmsM0))
    
    
    res$B1000 <- res$B/1000
    res$E1_10 <- res$E1/10
    res$E2_10 <- res$E2/10
    res$S1_10 <- res$S1/10
    res$S2_10 <- res$S2/10
    res$respO_10 <- res$respO/10
    cls <- c("B1000","E1","E2","respO","Mm","S1_10","S2_10")
    #cls <- c("B1000","E1_10","E2_10","respO_10","Mm","S1_100","S2_100")
    #cls <- c("B1000","E1_10","E2_10","respO","Mm","alpha")
    #cls <- c("B","E","respO","Mm","S1","S2")
    bo <- TRUE
    #bo <- res[,1] < 70
    matplot( res[bo,1], res[bo,cls], type="l", lty=1:10, col=1:10, ylim=c(0,70), xlab="time (d)", ylab="")
    legend("topright", inset=c(0.01,0.01), legend=cls, lty=1:10, col=1:10)
 
}

