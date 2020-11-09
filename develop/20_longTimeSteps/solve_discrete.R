adjust_parms_identity <- function(parms, x, dx, ts)
  list(is_adjusted = FALSE, parms = parms)
solveDiscreteAdj <- function(y, tout, func, parms, tstep
                             ,f_adjust_parms = adjust_parms_identity) {
  t <- seq(min(tout),max(tout),by = tstep)
  nstep <- length(t)
  x <- rbind(y, matrix(NA_real_, nrow = nstep, ncol = length(y)))
  outputs_list <- vector("list", nstep)
  for(i in 1:nstep) {
    resFunc <- func(t[i], x[i,], parms)
    ans_adjust <- f_adjust_parms(parms, x[i,], resFunc, tstep)
    if (isTRUE(ans_adjust$is_adjusted)){
      resFunc <- func(t[i], x[i,], ans_adjust$parms)
    }
    x[i+1,] <- x[i,] + resFunc[[1]] * tstep
    outputs_list[[i]] <- resFunc[[2]]
  }
  isOut <- c(t %in% tout)
  outputs_mat <- outputs_list[isOut] %>% bind_rows()
  ans <- cbind(time = t[isOut], x[c(isOut,FALSE),,drop = FALSE], outputs_mat)
  rownames(ans) <- NULL
  ans
}


solveDiscrete <- function(y, tout, func, parms, tstep) {
  t <- seq(min(tout),max(tout),by = tstep)
  nstep <- length(t)
  x <- rbind(y, matrix(NA_real_, nrow = nstep, ncol = length(y)))
  outputs_list <- vector("list", nstep)
  for(i in 1:nstep) {
    resFunc <- func(t[i], x[i,], parms)
    x[i+1,] <- x[i,] + resFunc[[1]] * tstep
    outputs_list[[i]] <- resFunc[[2]]
  }
  isOut <- c(t %in% tout)
  outputs_mat <- outputs_list[isOut] %>% bind_rows()
  ans <- cbind(time = t[isOut], x[c(isOut,FALSE),,drop = FALSE], outputs_mat)
  rownames(ans) <- NULL
  ans
}
