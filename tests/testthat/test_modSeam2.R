#require(testthat)
context("modSeam2")

parms0 <- list(
  cnB = 7.16
  ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
  #,cnE = 7.16
  ,cnIR = 4.5      ##<< between micr and enzyme signal
  ,cnIL = 30 ##<< N poor substrate
  #,kN = 0.05   ##<< (per day) enzyme turnover
  ,kN = 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
  ,kNB = 0.8     ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R)
  #,kR = 0.2      ##<< substrate decomposition rate N-rich (here simulate large N stock)
  #,kL = 1      ##<< substrate decomposition rate N-poor
  #,kR = 5e-3      ##<< substrate decomposition rate N-rich (here simulate large N stock)
  #,kL = 10e-3     ##<< substrate decomposition rate N-poor
  #,kL = 5e-2     ##<< substrate decomposition rate N-poor
  #,aE = 0.05   ##<< C-uptake allocated to enzymes
  #,kR = 1/(5*365)        ##<< 5 years
  #,kL = 1/(0.5*365)        ##<< 1/2 years
  ,kR = 1/(50)        ##<< 1/(x years)
  ,kL = 1/(1)        ##<< 1/(x years)
  #,aE = 0.003*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
  ,aE = 0.001*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
  ,km = 0.3     ##<< enzyme half-saturation constant
  #,km = 0.03     ##<< enzyme half-saturation constant
  #,km = 14     ##<< enzyme half-saturation constant
  ,m = 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
  ,tau = 1/60*365  ##<< biomass turnover rate (12 days)
  ,eps = 0.5      ##<< carbon use efficiency
  ,epsTvr = 0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
  ,iR = 0        ##<< input modelled explicitely
  ,iL = 300         # g/m2 input per year (half NPP)
  #,plantNUp = 300/70*1/4  # plant N uptake balancing N inputs
  ,plantNUp = 0	# organic N takeup by plants
  ,kIP = 10.57 #0.0289652*365          ##<< plant uptake iP I
  ,useFixedAlloc = FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
  ,iP = 10.57 #0.0289652*365          ##<< plant uptake iP I
  ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I
  ,iI = 0     ##<< input of mineral N
  ,l = 0.96   #0.00262647*365       ##<< leaching rate of mineralN l I
	,nu = 0.9     # microbial N use efficiency
)
parms0 <- within(parms0,{
            kmR <- kmL <- km
            eps1 <- eps2 <- eps
            cnER <- cnEL <- cnE
            kNR <- kNL <- kN
        })

parms <- parms0

x0 <- x0Orig <- c(  # aE = 0.003*365
        B = 1                     ##<< microbial biomass
        ,ER  = 2*parms0$km                  ##<< total enzyme pool
        ,EL  = 2*parms0$km                   ##<< total enzyme pool
        ,R = 1200                 ##<< N rich substrate
        ,RN = 1200/parms0$cnIR    ##<< N rich substrate N pool
        ,L = 170                  ##<< N poor substrate
        ,LN = 170/parms0$cnIL    ##<< N poor substrate N pool
        ,I =  0                    ##<< inorganic pool
)
x <- x0

x0 <- x0Orig <- c( #aE = 0.001*365
        B = 20                     ##<< microbial biomass
        ,ER  = 1.5*parms0$km                  ##<< total enzyme pool
        ,EL  = 1.5*parms0$km                   ##<< total enzyme pool
        ,R = 7000                 ##<< N rich substrate
        ,RN = 7000/parms0$cnIR    ##<< N rich substrate N pool
        ,L = 200                  ##<< N poor substrate
        ,LN = 200/parms0$cnIL    ##<< N poor substrate N pool
        ,I =  0                    ##<< inorganic pool
)
x <- x0


test_that("Regression to previous values", {
    #x0["L"] <- calcLSteady5( parms0$iL, x0["R"], parms=parms0 )
    #x0["B"] <- calcBSteady5( parms0$iL, x0["R"], L=x0["L"], parms=parms0 )
    #x0[c("ER","EL")] <- calcESteady5( x0["B"], parms=parms0 )

    parmsTvr0 <- within(parms0,{
                isTvrNil <- TRUE
                iR <- 160
            })
    parmsInit <- parms0
    #parmsInit <- within(parms0, {isFixedS <- TRUE})
    #parmsInit <- within(parms0, {isAlphaFix <- TRUE})
    #parmsInit <- within(parms0, {isAlphaMatch <- TRUE})
#    res <- derivSeam1(0, x0, parmsInit)
#    expect_equal( res, list(structure(c(-98.6526666666667, 1.80502623688156, 2.20997376311844,
#                                    -46.843, -13.3569658196672, 180, 6, 29.8401099297171), .Names = c("dB",
#                                    "dER", "dEL", "dR", "dRN", "dL", "dLN", "dI")), structure(c(0,
#                                    29.8401099297171, 17.9453240824173, 11.8947858472998, 0.472263868065967,
#                                    0.472263868065967, 0.823529411764706, 4.5, 30, 0.6, 0.6, 84,
#                                    120, 261.480666666667, 176.314, 85.1666666666667, 121.666666666667,
#                                    30.6849315068493, 43.8356164383562, 21.1385083713851, 4.5296803652968,
#                                    0.151905063590127, 303.005040860215, 46.028, 6.58305902624957,
#                                    6.42849162011173, 21.1595698924731), .Names = c("respO", "Mm",
#                                    "MmImb", "MmTvr", "alpha", "alphaC", "alphaN", "cnR", "cnL",
#                                    "limER", "limEL", "decR", "decL", "resp", "respB", "respTvr",
#                                    "tvrB", "revRC", "revLC", "revRN", "revLN", "pCsyn", "CsynReq",
#                                    "Csyn", "pNsyn", "NsynReq", "Nsyn")))
#    )
    times <- seq(0,2100, length.out = 101)
    #times <- seq(0,10000, length.out=101)
    res <- res1 <- as.data.frame(lsoda( x0, times, derivSeam2, parms = parmsInit))
    xE <- unlist(tail(res,1))
#    expect_equal( xE, xE0 <- structure(c(2100, 16.6605216752147, 0.393922360650999, 1.27212980956837,
#                            2785.23251444521, 408.593691017999, 370.747489822084, 12.3582496607361,
#                            22141.4895755274, 0, 10.0000681794826, 0.0914013079338929, 9.90866687154875,
#                            0.23644060760243, 0.215412880609298, 0.512430782984698, 6.81663122968414,
#                            30, 0.56767497776184, 0.809176062832645, 31.6221361139847, 299.99999411932,
#                            300.000460199995, 229.054405399706, 70.9460548002891, 101.351506857556,
#                            12.4849686746087, 52.2805297059652, 5.67779033179128, 5.40232140294974,
#                            0.99358436024401, 204.011880247985, 202.703013518371, 1.0064570659651,
#                            28.3104767483759, 14.2466396821218), .Names = c("time", "B",
##                            "ER", "EL", "R", "RN", "L", "LN", "I", "respO", "Mm", "MmImb",
#                            "MmTvr", "alpha", "alphaC", "alphaN", "cnR", "cnL", "limER",
#                            "limEL", "decR", "decL", "resp", "respB", "respTvr", "tvrB",
#                            "revRC", "revLC", "revRN", "revLN", "pCsyn", "CsynReq", "Csyn",
#                            "pNsyn", "NsynReq", "Nsyn"))
#    )
	expect_equal(xE, structure(c(2100, 14.2404542198333, 0.550230439121748, 0.873814982862216,
							2088.24808330046, 306.346054613259, 402.996631741352, 13.4332210580451,
							-938.842658269694, 43.5768876908186, 0, 1.53064233572481, 8.46935766565969,
							1.53064233572481, 10.0000000013845, 4.0166e-16, 0, 0.386385455567076,
							0.125169265206879, 0.386385455567076, 6.81663123077178, 30, 0.647154481660422,
							0.744423095308868, 27.0283821185336, 299.999999999945, 300.00000000926,
							239.359399123137, 60.6406008861234, 86.6294298373192, 13.4580635136764,
							94.0608509053946, 6.12032475866722, 9.7196212602241, 1.25151318537282,
							173.258859674592, 216.835747365411, 0.799032732285683, 30.2843222577389,
							12.0990823795106, -9.28955046219926e-09, 5.46833689440973e-11,
							-2.30926389122033e-11, -0.569999998615501, 331.186594750672,
							86.6294298372961, 331.186594750672, 15.3064233572481), .Names = c("time",
							"B", "ER", "EL", "R", "RN", "L", "LN", "I", "respO", "PhiB",
							"PhiU", "PhiTvr", "PhiBU", "PhiTotal", "immoPot", "MmImb", "alpha",
							"alphaC", "alphaN", "cnR", "cnL", "limER", "limEL", "decR", "decL",
							"resp", "respB", "respTvr", "tvrB", "revRC", "revLC", "revRN",
							"revLN", "pCsyn", "CsynReq", "Csyn", "pNsyn", "NsynReq", "Nsyn",
							"dR", "dL", "dB", "dI", "uC", "synB", "decC", "decN"))
							, tolerance = 1e-3
	)
    #res <- derivSeam1(0, xE, parmsInit)
    #res <- res1f <- as.data.frame(lsoda( x0, times, derivEezy5, parms=within(parms0, useFixedAlloc<-TRUE) ))

    #plotResSeam2(res, "topright", cls = c("B10","respO","PhiU","Rr","Lr","alpha100"))
})

