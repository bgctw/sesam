
deriv_inorg_nitrogen_ECA <-function(t, x, p){
  # plant uptake
  uINPlant_pot <- (p$kINPlant * x["IN"])/(TODO)
  uINPlant_req <- p$uINPlant # for simplicity
  plant_uptake <- min(uINPlant_pot, uINPlant_req)
  # microbial carbon part
  dec <- p$decS # for simplicity
  tvr <- p$tau * x["B"]
  synBC_C <- p$eps*(dec -tvr) # C available for biomass synthesis
  # microbial nitrogen part
  decN <- dec/p$cnIS
  tvrN <- tvr/p$cnB
  immoNPot <- p$iBN * x["IN"]
  synBC_N <- (decN + immoNPot - tvrN)*p$cnB # biomass synthesis constraint by N
  # growth limited either C or N
  synBC <- min(synBC_C, synBC_N) # biomass growth limited by either C or N
  nbal <- decN - synBC/p$cnB - tvrN # imbalance flux by mineral
  minN <- max(0, nbal) # if nbal is positive N is mineralized
  immoN <- max(0, -nbal) # if nbal is negative N is immobilized
  # inorganic N
  depoN <- p$iIN
  leachN <- p$lN * x["IN"]
  # here, we assume for simplicity that all tvr is mineralized
  minNPred <- tvrN
  dIN = depoN + minN +minNPred -leachN -immoN -plant_uptake
  dB = synBC
  # C mineralization
  respGr <- (1 - p$eps)/p$eps * synBC
  respTvr <- tvr
  respOverflow <- dec - dB - respGr - tvr
  resp <- respGr + respTvr + respOverflow
  # check mass balances
  meps <- .Machine$double.eps^(1/2)
  if(abs((dB)  - (p$decS - resp)) > meps)
    stop("mass balance C error")
  if(abs((dB/p$cnB + dIN)  -
         (p$decS/p$cnIS + depoN - leachN - plant_uptake)) > meps)
    stop("mass balance N error")
  list(
    dx = cbind(B = dB, IN = dIN)[1,],
    cbind(uINPlant_pot = uINPlant_pot, uINPlant_req = uINPlant_req,
          plant_uptake = plant_uptake,
          immoNPot = immoNPot, immoN = immoN,
          minN = minN, minNPred = minNPred,
          leachN = leachN, depoN = depoN,
          synBC = synBC, synBC_C = synBC_C, synBC_N = synBC_N)[1,]
  )
}

