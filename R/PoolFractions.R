
MultiPoolFractions_updateElement <- function(
  ### set states for given element
  .self       ##<< MultiPoolFractions object updated
  , x         ##<< numeric vector of state variables
  , element   ##<< scalar string of element to update
  , poolStatesMap ##<< named list with string vector of state variable names
    ## for each pool as returned by \code{\link{MultiPoolFractions_getPoolStatesMap}}
){
  ##seealso<< \code{\link{createMultiPoolFractions}},
  #message("MultiPoolFractions_updateElement for element ", element)
  poolNames <- names(poolStatesMap)
  fracEl <- lapply(poolStatesMap, function(extNames){
    structure(x[extNames], names = extNames )
  })
  .self$frac[poolNames] <- fracEl
  units <- .self$units[[element]]
  .self$poolPart[[element]] <- poolNames
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
  ##seealso<< \code{\link{createMultiPoolFractions}},
  ##value<< a list for each pool Name listing the state variable names
  poolStatesMap <- structure(lapply( poolNames, function(poolName){
    extNames <- paste(poolName, frac, sep = "_")
  }), names = poolNames)
}

MultiPoolFractions_updateScalars <- function(
  ### set states for given scalars, i.e. pools that are not devided into fractions
  .self       ##<< MultiPoolFractions object updated
  , x         ##<< numeric vector of state variables
  , poolNames ##<< string vector of pools that are scalars
){
  ##seealso<< \code{\link{createMultiPoolFractions}},
  # make sure list entry and scalar have the pool name
  .self$poolPart[["scalars"]] <- poolNames
  .self$frac[poolNames] <- structure(lapply(poolNames, function(poolName){
    structure(x[poolName], names = poolName)}), names = poolNames)
  .self$tot[poolNames] <- x[poolNames]
  .self$rel[poolNames] <- 1
  .self$units$scalars <- 1
  .self
}

sumMultiPoolFractions <- function(
  ### sum fraction-state variables across pools
  .self  ##<< MultiPoolFractions object mapping stateVars to pools
  , ds    ##<< data.frame with column names of state variables
  , keepVars = names(ds)    ##<< list of columns to keep
){
  resEl <- lapply( setdiff(names(.self$poolPart),"scalars"), function(element){
    units <- .self$units[[element]]
    pools <- .self$poolPart[[element]]
    unitM <- matrix(units, nrow = nrow(ds), ncol = length(units), byrow = TRUE)
    dsE <- as.data.frame(structure(do.call(cbind, lapply(pools, function(pool){
      rowSums( ds[, names(.self$frac[[pool]]), drop = FALSE]*unitM )
    })), dimnames = list(NULL, pools)))
  })
  keepVarsS <- union(keepVars, .self$poolPart$scalars)
  ##value<< a data.frame with columns named as in .self$tot
  cbind(ds[, keepVarsS, drop=FALSE], do.call(cbind, resEl))
}

### An (list) object as described by \code{\link{createMultiPoolFractions}}.
MulitPoolFractions <- list(
  ##describe<<
  className = "MultiPoolFractions"    ##<< string to identify the class
  , frac = list()   ##<< pool-list of state variable vector
    ##  with a named numeric vector of state variables
  , tot = numeric() ##<< named numeric vector with sums for each pool
  , rel = list()    ##<< pool-list of state variable vector normalized by total
  , stateVec = function(.self){ do.call(c, c(.self$frac, use.names = FALSE)) }
  #
  , units = list()
  , setX = function(.self, x) stop(
    "need to set specific function, see example createSesam4setX")
  , updateElement = MultiPoolFractions_updateElement  ##<< function
    ## to update frac and tot of fractions from given state vector
  , updateScalars = MultiPoolFractions_updateScalars
  ##end<<
)

createMultiPoolFractions <- function(
  ### construct a MultiPoolFractions object
  units   ##<< list for each element with a named numeric vector of fractions
  , setX  ##<< function to set state variable vector, see \code{\link{createSesam4setX}}
){
  ##details<< The example and its comments provide a good documentation on
  ## how to use a MultiPoolFractions object.
  pf <- MulitPoolFractions
  pf$units <- units
  pf$setX <- setX
  pf$poolParts <- list()  # mapping poolName -> partitioning, i.e. entry in units
  ##value<< a MultiPoolFractions object
  pf$setX(pf, numeric())
}
attr(createMultiPoolFractions,"ex") <- function(){
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppt
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
  )
  cnL = 30; cnB = 8; cnR = 10
  cpB = 47.3; cpR = 61; cpL = 919
  x0Vec <- c(B_C12 = 100, B_C13 = 27, B_C14 = 50
             , R_C12 = 10000, R_C13 = 2700, R_C14 = 5000
             , L_C12 = 1000, L_C13 = 270, L_C14 = 500
             , BN_N14 = 100/cnB, BN_N15 = 40/cnB
             , RN_N14 = 10000/cnR, RN_N15 = 600/cnR
             , LN_N14 = 1000/cnL, LN_N15 = 500/cnL
             , I_N14 = 10, I_N15 = 4
             , BP = 100/cpB, RP = 10000/cpR, LP = 1000/cpL, IP = 2
             , alpha = 0.5)
  #
  # create the Multipoolfractions object and set the state vector
  x <- createMultiPoolFractions(units, setX = createSesam4CNsetX(units))
  x <- x$setX(x, x0Vec)
  #
  # retrieve the original state vector back, names are sorted according to setX
  x$stateVec(x)
  #
  # totals across fractions for each element
  x$tot
  #
  # fractions for given pool (in units)
  x$frac[["B"]]
  # fractions in unit 1
  x$frac[["B"]] * x$units$C
  #
  # relative fractions to total (in units)
  # this is useful in mutiplying computed totals by its source
  x$rel[["B"]]
  # relative fractions sum to 1 at unit 1
  sum( x$rel[["B"]] * x$units$C )
}

createSesam4setX <- function(
  ### create a setter function for MultiPoolFfractions for Sesam4
  units  ##<< list with items C,N,P with a named numeric vector
    ## of fractions in these pools
){
  ##details<<
  ## returns a \code{function(.self,x) -> .self} to update a
  ## MulitPoolFractions object.
  ## It used \code{.self$updateElement}
  ## (\code{\link{MultiPoolFractions_updateElement}})
  ## and \code{.self$updateScalars}
  ## (\code{\link{MultiPoolFractions_updateScalars}})
  ## with arguments other than \code{.self} and \code{x}
  ## from local variables inside the create function (closure), that are
  ## created using function \code{\link{MultiPoolFractions_getPoolStatesMap}}.
  ##
  ## Look at the functions source code.
  #
  ##seealso<< \code{\link{createMultiPoolFractions}},
  ## \code{\link{createSesam4CNsetX}}
  .scalarPools <- c("alpha")
  required <- c("C","N","P")
  iMissing <- which( !(required %in% names(units)) )
  if (length(iMissing)) stop(
     "need to provide units for fractions in ", required[iMissing])
  .poolStatesMapC <- MultiPoolFractions_getPoolStatesMap(
    names(units$C), c("B","R","L"))
  .poolStatesMapN <- MultiPoolFractions_getPoolStatesMap(
    names(units$N), c("BN","RN","LN","I"))
  .poolStatesMapP <- MultiPoolFractions_getPoolStatesMap(
    names(units$P), c("BP","RP","LP","IP"))
  ##value<< a \code{function(.self,x) -> .self}
  ## that properly updates state of \code{.self}.
  function(.self,x){
    .self <- .self$updateElement(.self, x, "C", .poolStatesMapC)
    .self <- .self$updateElement(.self, x, "N", .poolStatesMapN)
    .self <- .self$updateElement(.self, x, "P", .poolStatesMapP)
    .self <- .self$updateScalars(.self, x, .scalarPools)
    .self
  }
}
attr(createSesam4setX,"ex") <- function(){
  # display source code
  createSesam4setX
}

createSesam4CNsetX <- function(
  ### create a setter function for MultiPoolFfractions for Sesam4 for CN isotopes
  units  ## list for C,N with a named numeric vector of fractions in these pools
){
  ##seealso<< \code{\link{createMultiPoolFractions}}
  .scalarPools <- c("BP","RP", "LP", "IP", "alpha")
  required <- c("C","N")
  iMissing <- which( !(required %in% names(units)) )
  if (length(iMissing)) stop(
    "need to provide units for fractions in ", required[iMissing])
  .poolStatesMapC <- MultiPoolFractions_getPoolStatesMap(
    names(units$C), c("B","R","L"))
  .poolStatesMapN <- MultiPoolFractions_getPoolStatesMap(
    names(units$N), c("BN","RN","LN","I"))
  ##value<< a function that properly updates frac and tot in .self
  function(.self,x){
    .self <- .self$updateElement(.self, x, "C", .poolStatesMapC)
    .self <- .self$updateElement(.self, x, "N", .poolStatesMapN)
    .self <- .self$updateScalars(.self, x, .scalarPools)
    .self
  }
}






