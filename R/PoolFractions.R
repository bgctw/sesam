
MultiPoolFractions_updateElement <- function(
  ### set states for given element
  .self       ##<< MultiPoolFractions object updated
  , x         ##<< numeric vector of state variables
  , element   ##<< scalar string of element to update
  , poolStatesMap ##<< named list with string vector of state variable names
    ## for each pool, see \code{\link{MultiPoolFractions_getPoolStatesMap}}
){
  #message("MultiPoolFractions_updateElement for element ", element)
  poolNames <- names(poolStatesMap)
  fracEl <- lapply(poolStatesMap, function(extNames){
    structure(x[extNames], names = extNames )
  })
  .self$frac[poolNames] <- fracEl
  units <- .self$units[[element]]
  .self$tot[poolNames] <- sapply(poolNames, function(poolName){
    sum(fracEl[[poolName]]*units)
  })
  .self$rel[poolNames] <- lapply(poolNames, function(poolName){
    fracEl[[poolName]]/.self$tot[poolName]
  })
  .self
}

MultiPoolFractions_getPoolStatesMap <- function(
  ### construct state variable names for pools devided into fractions
  frac        ##<< numeric vector of fractions for element pools
  , poolNames ##<< numeric vector of pools that are devided into these fractions
){
  ##value<< a list for each pool Name listing the state variable names
  poolStatesMap <- structure(lapply( poolNames, function(poolName){
    extNames <- paste(poolName, frac, sep = "_")
  }), names = poolNames)
}



MultiPoolFractions_updateScalars <- function(
  ### set states for given scalars, i.e. pools that are not devided into fractions
  .self       ##<< MultiPoolFractions object updated
  , x         ##<< numeric vector of state variables
  , poolNames ##<< numeric vector of pools that are scalars
){
  # make sure list entry and scalar have the pool name
  .self$frac[poolNames] <- structure(lapply(poolNames, function(poolName){
    structure(x[poolName], names = poolName)}), names = poolNames)
  .self$tot[poolNames] <- x[poolNames]
  .self$rel[poolNames] <- 1
  .self
}

### An object (list) with list \code{frac} and vector \code{tot}
MulitPoolFractions <- list(
  ##describe<<
  className = "MultiPoolFractions"    ##<< string to identify the class
  , frac = list()   ##<< pool-list of state variable vector
    ##  with a named numeric vector of state variables
  , tot = numeric() ##<< named numeric vector with sums for each pool
  , rel = list()    ##<< pool-list of state variable vector normalized by total
  , updateElement = MultiPoolFractions_updateElement  ##<< function
    ## to update frac and tot of fractions from given state vector
  , updateScalars = MultiPoolFractions_updateScalars
  , setX = function(.self, x) stop(
    "need to set specific function, see createSesam4setX")
  , poolNames = function(.self){ names(.self$frac) }
  , stateNames = function(.self){
    do.call(c, sapply(.self$frac, function(fraci) names(fraci)))  }
  ##end<<
)

createMultiPoolFractions <- function(
  ### construct a MultiPoolFractions object
  units   ## list for each element with a named numeric vector of fractions
  , setX  ##<< function to set state variable vector, see \code{\link{createSesam4setX}}
){
  pf <- MulitPoolFractions
  pf$units <- units
  pf$setX <- setX
  pf
}



createSesam4setX <- function(
  ### create a setter function for MultiPoolFfractions for Sesam4
  units  ## list for C,N,P with a named numeric vector of fractions in these pools
){
  .scalarPools <- c("alpha")
  required <- c("C","N","P")
  iMissing <- which( !(required %in% names(units)) )
  if (length(iMissing)) stop(
     "need to provide units for fractions in ", required[iMissing])
  .poolStatesMapC <- MultiPoolFractions_getPoolStatesMap(
    names(units$C), c("B","R","L"))
  .poolStatesMapN <- MultiPoolFractions_getPoolStatesMap(
    names(units$N), c("RN","LN","I"))
  .poolStatesMapP <- MultiPoolFractions_getPoolStatesMap(
    names(units$P), c("RP","LP","IP"))
  ##value<< a function that properly updates frac and tot in .self
  function(.self,x){
    .self <- .self$updateElement(.self, x, "C", .poolStatesMapC)
    .self <- .self$updateElement(.self, x, "N", .poolStatesMapN)
    .self <- .self$updateElement(.self, x, "P", .poolStatesMapP)
    .self <- .self$updateScalars(.self, x, .scalarPools)
    .self
  }
}

createSesam4CNsetX <- function(
  ### create a setter function for MultiPoolFfractions for Sesam4 for CN isotopes
  units  ## list for C,N with a named numeric vector of fractions in these pools
){
  .scalarPools <- c("RP", "LP", "IP", "alpha")
  required <- c("C","N")
  iMissing <- which( !(required %in% names(units)) )
  if (length(iMissing)) stop(
    "need to provide units for fractions in ", required[iMissing])
  .poolStatesMapC <- MultiPoolFractions_getPoolStatesMap(
    names(units$C), c("B","R","L"))
  .poolStatesMapN <- MultiPoolFractions_getPoolStatesMap(
    names(units$N), c("RN","LN","I"))
  ##value<< a function that properly updates frac and tot in .self
  function(.self,x){
    .self <- .self$updateElement(.self, x, "C", .poolStatesMapC)
    .self <- .self$updateElement(.self, x, "N", .poolStatesMapN)
    .self <- .self$updateScalars(.self, x, .scalarPools)
    .self
  }
}






