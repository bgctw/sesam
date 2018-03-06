from sympy import Eq, Symbol, symbols, solve, simplify, assuming, Q, refine

dL, dR, B, alpha, aE, kN, km, tau, eps, synE = symbols('dL, dR, B, alpha, aE, kN, km, tau, eps, synE')
m = symbols('m')
betaL, betaR, betaB, uImmPot = symbols('betaL, betaR, betaB, uImmPot')

EL = (1-alpha)*aE*B/kN
decL = dL * (EL/(km + EL)).simplify()

ER = (alpha)*aE*B/kN
decR = dR * (ER/(km + ER)).simplify()

uC = (decL + decR).simplify()

synBC = eps*(uC - m*B).simplify()

tvrB = tau*B

eqB = Eq(synBC, tvrB)

'with assuming(alpha >= 0, alpha <= 1, B > 0, aE > 0, kN > 0, km > 0, eps > 0, m >= 0, tau > 0):
ans1 = solve(eqB, B)
ans2 = refine(ans1, Q.positive(B))
Bs = ans2[2].simplify()
print(Bs)


uN = (decL/betaL + decR/betaR + uImmPot).simplify()
synBN = (betaB*uN/eps).simplify()

eqBN = Eq(synBN, tvrB)

# takes long
#with assuming(alpha >= 0, alpha <= 1, B > 0, aE > 0, kN > 0, km > 0, eps > 0, m >= 0, tau > 0):
ans1 = solve(eqBN, B)
ans2 = refine(ans1N, Q.positive(B))
BsN = ans2[2].simplify()
print(Bs)
