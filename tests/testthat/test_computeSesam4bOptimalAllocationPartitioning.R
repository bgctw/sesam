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
