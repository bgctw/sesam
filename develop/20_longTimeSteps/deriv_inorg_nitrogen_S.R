deriv_inorg_nitrogen_S <-function(t, x, p){
  # plant uptake
  uINPlant_pot <- p$kINPlant * x["IN"]
  uINPlant_req <- p$uINPlant # for simplicity
  plant_uptake <- min(uINPlant_pot, uINPlant_req)
  # microbial carbon part
  dec <- p$kS * x["SC"]
  tvr <- p$tau * x["B"]
  synBC_C <- p$eps*(dec -tvr)
  # microbial nitrogen part
  cnS <- x["SC"]/x["SN"]
  decN <- dec/cnS
  tvrN <- tvr/p$cnB
  immoNPot <- p$iB * x["IN"]
  synBC_N <- (decN + immoNPot - tvrN)*p$cnB
  # growth limited either C or N
  synBC <- min(synBC_C, synBC_N) # biomass growth limited by either C or N
  nbal <- decN - synBC/p$cnB - tvrN
  minN <- max(0, nbal) # if nbal is positive N is mineralized
  immoN <- max(0, -nbal) # if nbal is negative N is immobilized
  # inorganic N
  depoN <- p$iIN
  leachN <- p$lN * x["IN"]
  # a part of turnover is respired and a corresponding part of N is mineralized
  minNPred <- p$epsTvr * tvrN
  dB = synBC - tvr
  dSC = +p$iS - dec + p$epsTvr*tvr
  dSN = +p$iS/p$cnIS - decN  + p$epsTvr*tvrN
  dIN = depoN + minN +minNPred -leachN -immoN
  respG = (1 - p$eps)*synBC
  respTvr = (1 - p$epsTvr)*tvr
  # check mass balances
  meps <- .Machine$double.eps^(1/2)
  if(abs((dSC + dB)  - (p$iS - respG - respTvr)) > meps)
    stop("mass balance C error")
  if(abs((dSN + dB/p$cnB)  - (p$iS/p$cnIS + depoN - leachN -plant_uptake)) > meps)
    stop("mass balance N error")
  list(
    dx = cbind(B = dB, SC = dSC, SN = dSN, IN = dIN)[1,],
    cbind(uINPlant_pot = uINPlant_pot, uINPlant_req = uINPlant_req,
          immoNPot = immoNPot, immoN = immoN,
          minN = minN, minNPred = minNPred,
          leachN = leachN, depoN = depoN,
          synBC = synBC, synBC_C = synBC_C, synBC_N = synBC_N)[1,]
  )
}
