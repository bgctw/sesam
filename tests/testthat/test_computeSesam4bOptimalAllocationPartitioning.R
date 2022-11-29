#require(testthat)
context("computeSesam4bOptimalAllocationPartitioning")

test_that("computeSesam4bOptimalAllocationPartitioning", {
  B=1
  synB = 15
  params = within(list(
    aE = 0.1,
    e_P = 0,
    tau = 365/30
  ), kmN <- aE*B/2)
  nstep = 200
  dL = 0.7
  dR = 0.5
  dPs = seq(0,1, length.out=nstep)
  alpha3 = sapply(dPs, function(dP){computeSesam4bOptimalAllocationPartitioning(dL,dR,dP,params, B, synB)} )
  .tmp.f <- function(){
    plot(dPs, alpha3["P",], type="l", ylim=c(-0.1,1.1), ylab=expression(alpha), xlab=expression(d[P]))
    #plot(dPs, alpha3["R",], type="l",  ylab=expression(alpha), xlab=expression(d[P]))
    lines(dPs, alpha3["L",], type="l", lty="dashed")
    lines(dPs, alpha3["R",], type="l", lty="dotted")
  }
  expect_true(all.equal(colSums(alpha3), rep(1,nstep)))
  expect_true( all(alpha3["L",] > alpha3["R",]))
  expect_true( alpha3["P",1] == 0) # no enzyme allocation to P for dP=0
  expect_true( alpha3["P",nstep] > alpha3["L",nstep])
  #
  p = within(params, B<-B, synB <- synB)
  uL <- u_decomp(dL, alpha3["L",], p)
  uR <- u_decomp(dR, alpha3["R",], p)
  uP <- u_biomin(dPs, alpha3["P",], p)
  # not equal, but derivatives need to be equal expect_equal(uL, uR)
  du <- du_dalpha(t(alpha3), dL, dR, dPs, p)
  expect_equal(du[,"L"], du[,"R"])
  is_Pallocated <- alpha3["P",] > 0
  expect_equal(du[is_Pallocated,"L"], du[is_Pallocated,"P"])
})

.tmp.check_time_evolution <- function(){
  B <- 1
  parms = within(list(
    aE = 0.1,
    e_P = 0,
    tau = 365/30
    ,cnIR = 4.5     ##<< between micr and enzyme signal
    ,cnIL = 30      ##<< N poor substrate
    , cpE = 50
    , cpB = 40
    #, cpBW = 50
    , cpIR = 40
    , cpIL = 40*3
  ), kmN <- aE*B/2)

  #library(deSolve)
  alpha0 = alpha <- c(L=1,R=1,P=1)/3 #compute_alpha3_relative(dL,dR,dP,p)
  dLPot <- dLPPot <- 1.0
  dRPot <- dRPPot <- 0.7
  synB <- parms$tau * B # steady state B
  limE <- c(C=0.8, N=0.1, P=0.1)
  cnL <- parms$cnIL
  cnR <- parms$cnIR
  cpL <- parms$cpIL
  cpR <- parms$cpIR

  da0 <- calc_dAlphaP_propto_du(
    alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
    cnL, cnR, cpL, cpR
  )
  da0_opt <- calc_dAlphaP_optimal(
    alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
    cnL, cnR, cpL, cpR
  )
  (da0_opt$alphaTarget - alpha)
  da0$dud

  times <- seq(0,0.2,length.out=201)
  #times <- seq(0,2,length.out=201)
  #times <- seq(0,0.03,length.out=201)
  res_ode <- lsode(alpha0, times, function(t,alpha,parms){
    calc_dAlphaP_propto_du(
      alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
      cnL, cnR, cpL, cpR
    )
    }, parms)
  res_ode_opt <- lsode(alpha0, times, function(t,alpha,parms){
    calc_dAlphaP_optimal(
      alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
      cnL, cnR, cpL, cpR
    )
  }, parms)
  plot(times, res_ode[,"L"], type="l")
  lines(times, res_ode_opt[,"L"], col="blue")
  # currently the derivative-based dynamics is much faster than the optimal

  plot(times, res_ode[,"P"], type="l")
  lines(times, res_ode_opt[,"P"], col="blue")

  plot(times, res_ode[,"R"], type="l")
  lines(times, res_ode_opt[,"R"], col="blue")


  alpha <- res_ode_opt[length(times),2:4]

}
