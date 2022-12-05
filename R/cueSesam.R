cueSesam <- function(
  ### compute carbon use efficiency from sesam results
  resSesam    ##<< data.frame of seasam output columns
){
  # ##value<< input data.frame with computed columns
  # ans <- resSesam %>% mutate(
  #   cueDef = synB / (uptakeC)  ##<<
  #     ## biomass synthesis / uptake
  #   , cueDB = (synB - tvrB - tvrBPred) / (synB - tvrB - tvrBPred + resp) ##<<
  #     ## biomass change / (biomass change + resp )
  #   , cueSyn = (synB ) / (synB + resp) ##<<
  #     ## biomass synthesis / (biomass synthesis + resp )
  # )
  # ans
  ##value<< input data.frame with computed columns
  a <- resSesam
  # different variants of SESAM have different returns
  if (is_null(a$respTotal)) a$respTotal = a$resp
  if (is_null(a$tvrBPred)) a$tvrBPred = 0.0
  synB0 <- pmax(0, a$synB)  # non-negative biomass synthesis
  a$cueDef =  synB0 / (a$uptakeC)  ##<<
  ## biomass synthesis / uptake
  a$cueDB = (a$synB - a$tvrB - a$tvrBPred) / (a$synB - a$tvrB - a$tvrBPred + a$respTotal) ##<<
  # ## biomass change / (biomass change + resp )
  a$cueSyn = (synB0 ) / (synB0 + a$respTotal) ##<<
  ## biomass synthesis / (biomass synthesis + resp )
  a
}
