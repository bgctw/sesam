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
  a$cueDef =  a$synB / (a$uptakeC)  ##<<
  ## biomass synthesis / uptake
  a$cueDB = (a$synB - a$tvrB - a$tvrBPred) / (a$synB - a$tvrB - a$tvrBPred + a$resp) ##<<
  ## biomass change / (biomass change + resp )
  a$cueSyn = (a$synB ) / (a$synB + a$resp) ##<<
  ## biomass synthesis / (biomass synthesis + resp )
  a
}
