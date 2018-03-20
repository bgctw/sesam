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
  synB0 <- pmax(0, a$synB)  # non-negative biomass synthesis
  a$cueDef =  synB0 / (a$uptakeC)  ##<<
  ## biomass synthesis / uptake
  a$cueDB = (a$synB - a$tvrB - a$tvrBPred) / (a$synB - a$tvrB - a$tvrBPred + a$resp) ##<<
  # ## biomass change / (biomass change + resp )
  a$cueSyn = (synB0 ) / (synB0 + a$resp) ##<<
  ## biomass synthesis / (biomass synthesis + resp )
  a
}