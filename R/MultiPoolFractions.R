createMultiPoolFractions <- function(
  ### construct a MultiPoolFractions object
  units   ##<< list for each element with a named numeric vector of fractions
  , setX  ##<< function to set state variable vector, see \code{\link{createSesam4setX}}
){
  ##details<< The example and its comments provide a good documentation on
  ## how to use a MultiPoolFractions object.
  #
  ##seealso<<
  ## \itemize{
  ## \item sum pools of state variables: \code{\link{sumMultiPoolFractions}}
  ## \item sum pools of output fractions: \code{\link{sumMultiPoolFractionsVars}}
  ## \item set all state variables of one pool: \code{\link{setMultiPoolFractionsPool}}
  ## \item set all state variables of elemental pools of one entity:
  ##    \code{\link{setMultiPoolFractionsElements}}
  ## }
  ##value<< a MultiPoolFractions object with items
  pf <- list(
    frac = list()      ##<< pool-list of state variable vector
    ##  with a named numeric vector of state variables
    , tot = numeric()   ##<< named numeric vector with sums for each pool
    , rel = list()      ##<< pool-list of state variable vector normalized by total
    , poolPart = list() ##<< mapping element -> pool
    ## Which partitioning, i.e. element, is applied to which pools.
    ## Also contains entry "scalars" listing pools with a sole single fraction.
    , stateVec = function(.self){ do.call(c, c(.self$frac, use.names = FALSE)) }
    , units = units     ##<< from argument list
    , setX = setX       ##<< from argument list
    , updateElement = MultiPoolFractions_updateElement  ##<< function
    ## to update frac and tot of fractions from given state vector
    , updateScalars = MultiPoolFractions_updateScalars
    ##end<<
  )
  class(pf) <- c("MultiPoolFractions","list")
  pf$setX(pf, numeric())
}
attr(createMultiPoolFractions,"ex") <- function(){
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppt
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
    , P = c(tot = 1)
  )
  x <- createMultiPoolFractions(units, setX = createSesam4setX(units))
  #
  # initializing state variables by elemental ratios and atomic ratios
  cnB = 8;  cnR = 10; cnL = 30
  cpB = 47.3; cpR = 61; cpL = 919
  aR <- list(
    C = c(C13 = 27, C14 = 120) # ~ C3 plants and atmosphere 1990
    , N = c(N14 = 0.4))
  x0Vec <- x$stateVec(x)
  x0Vec["alpha"] <- 0.5    # scalars
  x0Vec["IP_tot"] <- 2
  x0Vec <- setMultiPoolFractionsPool(  # single element pools
    x, x0Vec, "I", 10, rel = aR$N, element = "N")
  x0Vec <- setMultiPoolFractionsElements( # multiple element pools
    x, x0Vec, "B", 100, ce = c(cnB, cpB), rel = aR)
  x0Vec <- setMultiPoolFractionsElements(
    x, x0Vec, "R", 10000, ce = c(cnR, cpR), rel = aR)
  x0Vec <- setMultiPoolFractionsElements(
    x, x0Vec, "L", 1000, ce = c(cnL, cpL), rel = aR)
  x <- x$setX(x, x0Vec)  # setting the state vector will compute totals and rel
  #
  # retrieve the original state vector back, names are sorted according to setX
  x$stateVec(x)
  #
  # totals across fractions for each element
  x$tot
  x$tot["BC"]/x$tot["BN"]  # equals cnB
  #
  # fractions for given pool (in units)
  x$frac[["BC"]]
  # fractions in unit 1
  x$frac[["BC"]] * x$units$C
  #
  # relative fractions to total (in units)
  # this is useful in mutiplying computed totals by its source
  x$rel[["BC"]]
  # relative fractions sum to 1 at unit 1
  sum( x$rel[["BC"]] * x$units$C )
  # atomic ratios
  x$rel[["BC"]]  / x$rel[["BC"]][1]
}

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
  #.self$units$scalars <- 1 # better do not modify units
  .self
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
    names(units$C), c("BC","RC","LC","resp"))
  .poolStatesMapN <- MultiPoolFractions_getPoolStatesMap(
    names(units$N), c("BN","RN","LN","I","leachN"))
  .poolStatesMapP <- MultiPoolFractions_getPoolStatesMap(
    names(units$P), c("BP","RP","LP","IP","leachP"))
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
  .scalarPools <- c("BP","RP", "LP", "IP", "leachP", "alpha")
  required <- c("C","N")
  iMissing <- which( !(required %in% names(units)) )
  if (length(iMissing)) stop(
    "need to provide units for fractions in ", required[iMissing])
  .poolStatesMapC <- MultiPoolFractions_getPoolStatesMap(
    names(units$C), c("BC","RC","LC","resp"))
  .poolStatesMapN <- MultiPoolFractions_getPoolStatesMap(
    names(units$N), c("BN","RN","LN","I","leachN"))
  ##value<< a function that properly updates frac and tot in .self
  function(.self,x){
    .self <- .self$updateElement(.self, x, "C", .poolStatesMapC)
    .self <- .self$updateElement(.self, x, "N", .poolStatesMapN)
    .self <- .self$updateScalars(.self, x, .scalarPools)
    .self
  }
}


sumMultiPoolFractions <- function(
  ### sum fraction-state variables across pools
  .self  ##<< MultiPoolFractions object mapping stateVars to pools
  , ds    ##<< data.frame with column names of state variables
  , keepVars = names(ds)    ##<< list of columns to keep
){
  if (!(is.matrix(ds) || is.data.frame(ds))) ds <- as.data.frame(t(ds))
  ##seealso<< \code{\link{createMultiPoolFractions}}
  resEl <- lapply( setdiff(names(.self$poolPart),"scalars"), function(element){
    units <- .self$units[[element]]
    pools <- .self$poolPart[[element]]
    unitM <- matrix(units, nrow = nrow(ds), ncol = length(units), byrow = TRUE)
    dsE <- as.data.frame(structure(do.call(cbind, lapply(pools, function(pool){
      rowSums( ds[, names(.self$frac[[pool]]), drop = FALSE]*unitM )
    })), dimnames = list(NULL, pools)))
  })
  sumDs <- do.call(cbind, resEl)
  # at least keep scalars, but at least drop columns that will be added by sumDS
  keepVarsS <- setdiff( union(keepVars, .self$poolPart$scalars), names(sumDs))
  ##value<< a data.frame with columns named as in .self$tot
  cbind(ds[, keepVarsS, drop = FALSE], sumDs)
}


sumMultiPoolFractionsVars <- function(
  ### sum fraction variables across pools
  .self  ##<< MultiPoolFractions object mapping stateVars to pools
  , ds    ##<< data.frame with column names of state variables
  , vars  ##<< list of element (corresponding to item in \code{.self$units})
    ## -> string vector of variables
  , keepVars = names(ds)    ##<< list of columns to keep
){
  ##seealso<< \code{\link{createMultiPoolFractions}}
  if (!(is.matrix(ds) || is.data.frame(ds))) ds <- as.data.frame(t(ds))
  resEl <- lapply( names(vars), function(element){
    units <- .self$units[[element]]
    pools <- vars[[element]] #outer(vars[[element]], names(units), paste, sep = "_")
    unitM <- matrix(units, nrow = nrow(ds), ncol = length(units), byrow = TRUE)
    dsE <- as.data.frame(structure(do.call(cbind, lapply(pools, function(pool){
      fracNames <- paste(pool, names(units), sep = "_")
      rowSums( ds[, fracNames, drop = FALSE]*unitM )
    })), dimnames = list(NULL, pools)))
  })
  ##value<< a data.frame with columns of sums
  cbind(ds[, keepVars, drop = FALSE], do.call(cbind, resEl))
}
attr(sumMultiPoolFractionsVars,"ex") <- function(){
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppTrillion
    , N = c(N14 = 1, N15 = 0.01)) # 15N in percent
  x <- createMultiPoolFractions(units, setX = createSesam4CNsetX(units))
  ds <- data.frame(time = 1:5, decNL_N14 = 1:5, decNL_N15 = 6:10)
  dsSum <- sumMultiPoolFractionsVars(x, ds, vars = list(N = c("decNL")))
  dsSum
}

setMultiPoolFractionsPool <- function(
  ### set state variables corresponding to pool by fractions
  .self       ##<< MultiPoolFractions object mapping stateVars to pools
  , xvec      ##<< named numeric vector of state variables
  , pool      ##<< scalar string: pool to set
  , firstFrac ##<< scalar numeric: value of first fraction
  , total = NULL ##<< alternative to firstFrac: value of the sum of all fractions.
    ## If specified this will compute and overide firstFrac.
  , rel = numeric(0)    ##<< a named vector giving relative
    ## amount of other fractions to the first fraction (usually the heavy isotope).
    ## If only one vector is given, it is repliaced for all elements.
    ## If not specified or zero length, then the current proportions used.
  , element = substr(pool, nchar(pool), nchar(pool)) ##<< the element that
    ## for wich this pool tracks fractions
){
  ##seealso<< \code{\link{createMultiPoolFractions}},
  stateNames <- paste(pool, names(.self$units[[element]]), sep = "_")
  # if missing rel, then use the current partitioning
  if (!length(rel)) {
      xp <- xvec[stateNames]
      rel <- if (length(xp) == 1L) c() else xp[-1]/xp[1]
  }
  if (length(total)) {
    firstFrac = total/sum(c(1,rel)*.self$units[[element]])
  }
  if (length(rel) == 0) {
    # single fraction
    xvec[stateNames] <- firstFrac
  } else {
    relA <- c(1, rel)
    xvec[stateNames] <- firstFrac*relA
  }
  ##value<< state vector xvec with all entries for given pool updated
  xvec
}

setMultiPoolFractionsElements <- function(
  ### set state variables corresponding to pool by fractions
  .self   ##<< MultiPoolFractions object mapping stateVars to pools
  , xvec  ##<< named numeric vector of state variables
  , poolBase  ##<< scalar string: part of the poolname without element
  , value ##<< scalar numeric: value of first fraction of the first element (C)
  , ce = numeric()   ##<< named numeric vector of C:Element ratios, of others elements
  ## If not specified or zero length, then the current elemental ratios are used.
  , rel = list()    ##<< list: for each element a named vector giving relative
  ## amount of other fractions to the first fraction (usually the heavy isotope).
  ## If only one vector is given, it is repliaced for all elements.
  ## If not specified or zero length, then the current proportions used.
){
  ##seealso<< \code{\link{createMultiPoolFractions}},
  ## \code{\link{setMultiPoolFractionsPool}}
  elements <- names(.self$units)
  # if missing ce, then use the current elemental partitioning
  if (!length(ce) && length(elements != 1L)) {
    xs <- sumMultiPoolFractions(.self, xvec)
    poolNames <- paste0(poolBase,elements)
    xsp <- unlist(xs[poolNames])
    ce <- setNames(xsp[1L]/xsp[-1L], elements[-1])
  }
  # if only one duplicate vector to all list entries if
  if (!is.list(rel)) rel <- setNames(
    lapply(elements, function(element) rel), elements)
  iE <- 1L
  element <- elements[iE]
  pool <- paste0(poolBase, element)
  xvec <- setMultiPoolFractionsPool(
    .self, xvec, pool, firstFrac = value, rel = rel[[element]], element = element)
  stateNames <- paste(pool, names(.self$units[[element]]), sep = "_")
  cTot <- sum(xvec[stateNames]*.self$units[[iE]])
  eTot <- cTot / c(1, ce)
  while (iE < length(elements)) {
    iE <- iE + 1L
    element <- elements[iE]
    pool <- paste0(poolBase, element)
    xvec <- setMultiPoolFractionsPool(
      .self, xvec, pool, firstFrac = NULL
      , total = eTot[iE], rel = rel[[element]], element = element)
  }
  ##value<< state vector xvec with all entries for pools updated
  xvec
}

elementMultiPoolFractions <- function(
  ### get a matrix of state-variable of pools of one partitioning
  .self      ##<< MultiPoolFractions object
  , element  ##<< string scalar of element/partitioning
  , pools = .self$poolPart[[element]]
){
  ans <- do.call(cbind, sapply(
    pools, function(pool) .self$frac[[pool]], simplify = FALSE))
  rownames(ans) <- names(.self$units[[element]])
  ##value<< named numeric matrix with pools in columns and fraction in rows
  ans
}








