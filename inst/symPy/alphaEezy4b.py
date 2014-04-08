from __future__ import division 
from sympy import symbols, solve, Eq, diff, sqrt, simplify 

#E1, E2 = symbols('E1 E2',real=True, positiv=True)
alpha = symbols('alpha' ,real=True, positiv=True)
d1p, d2p = symbols('d1p d2p',real=True, positiv=True)       # potential decomposition fluxes kS 
kN1, kN2, km1, km2 = symbols('kN1 kN2 km1 km2',real=True, positiv=True)
aE, B = symbols('aE B',real=True, positiv=True)
cnS1, cnS2, cnE, cnB = symbols('cnS1 cnS2 cnE cnB',real=True, positiv=True)

E1 = alpha * aE * B / kN1
E2 = (1-alpha) * aE * B / kN2

#------------ alpha C-limited
e1 = d1p/(kN1*(E1+km1))      # efficiencies for C acquiring: decomp/E_tvr, E1 cancels
e2 = d2p/(kN2*(E2+km2))
eq1 = Eq( alpha - e1/(e1+e2))
#res = solve(eq1, alpha)
#alphaC = res[0]
alphaC = (2*B*aE*d1p + d1p*kN2*km2 + d2p*kN1*km1 - sqrt(4*B**2*aE**2*d1p*d2p + 4*B*aE*d1p*d2p*kN1*km1 + 4*B*aE*d1p*d2p*kN2*km2 + d1p**2*kN2**2*km2**2 + 2*d1p*d2p*kN1*kN2*km1*km2 + d2p**2*kN1**2*km1**2))/(2*B*aE*d1p - 2*B*aE*d2p)

#----------- alpha N-limited
e1N = d1p/(kN1*(E1+km1)) *  cnE/cnS1   
e2N = d2p/(kN2*(E2+km2)) *  cnE/cnS2    
eq1N = Eq( alpha - e1N/(e1N+e2N))
#resN = solve(eq1N, alpha)
#alphaN = resN[0]
alphaN = -(2*B*aE*cnS2*d1p + cnS1*d2p*kN1*km1 + cnS2*d1p*kN2*km2 - sqrt(4*B**2*aE**2*cnS1*cnS2*d1p*d2p + 4*B*aE*cnS1*cnS2*d1p*d2p*kN1*km1 + 4*B*aE*cnS1*cnS2*d1p*d2p*kN2*km2 + cnS1**2*d2p**2*kN1**2*km1**2 + 2*cnS1*cnS2*d1p*d2p*kN1*kN2*km1*km2 + cnS2**2*d1p**2*kN2**2*km2**2))/(2*B*aE*cnS1*d2p - 2*B*aE*cnS2*d1p)




#------------ biomass C limited
m, tau, eps = symbols('m tau eps',real=True, positiv=True)   # maintenance and microbial turnover
d1 = d1p*E1/(km1+E1)
d2 = d2p*E2/(km2+E2)
eqB = Eq( d1 + d2 - aE/eps*B - m*B - tau/eps*B )  
eqBC = eqB.subs( alpha, alphaC)
#resBC = solve(eqBC, B)    # takes long time
#--- N limited
eqBN0 = Eq( d1/cnS1 + d2/cnS2 - aE*B/cnE - tau*B/cnB )  # N balance without respiration 
eqBN = eqBN0.subs( alpha, alphaN)
#resBN = solve(eqBN, B)    # takes long time

#------------- optimize aE for maximum growth
fGrowth = d1 + d2 - aE/eps*B - m*B - tau/eps*B
dFGrowth = diff( fGrowth, aE)
eqdG = Eq(dFGrowth)
eqdGC = eqdG.subs( alpha, alphaC)
#res = solve(eqdGC, aE)
#aEC = -(kN1*km1 + kN2*km2)/B
eqdGN = eqdG.subs( alpha, alphaN)
#res = solve(eqdGN, aE)

print(fGrowth)





