from __future__ import division 
from sympy import symbols, solve, Eq, diff, sqrt, simplify 

#E1, E2 = symbols('E1 E2',real=True, positiv=True)
alpha = symbols('alpha' ,real=True, positiv=True)
S1, S2, kS1, kS2 = symbols('S1 S2 kS1 kS2',real=True, positiv=True)       # potential decomposition fluxes kS
kN1, kN2, km1, km2 = symbols('kN1 kN2 km1 km2',real=True, positiv=True)
aE,B = symbols('aE B',real=True, positiv=True)
cnS1, cnS2, cnE, cnB = symbols('cnS1 cnS2 cnE cnB',real=True, positiv=True)
i2 = symbols('i2',real=True, positiv=True)

E1 = alpha * aE * B / kN1
E2 = (1-alpha) * aE * B / kN2


#------------ biomass C limited
m, tau, eps = symbols('m tau eps',real=True, positiv=True)   # maintenance and microbial turnover
d1p = kS1*S1
d2p = kS2*S2
d1 = d1p*S1*E1/(km1+E1)
#d2 = kS2*S2*E2/(km2+E2)
d2 = i2
eqB = Eq( d1 + d2 - aE/eps*B - m*B - tau/eps*B )  
resB = solve(eqB, B)    
Beq = resB[1]

#------------- 
E2B = E1.subs(B, Beq)
res = solve( Eq(kS2*S2*E2B/(km2+E2B) -i2), S2 )
S2eq = res[0].simplify()

