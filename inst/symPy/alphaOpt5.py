from __future__ import division 
from sympy import symbols, solve, Eq, diff, sqrt, simplify 

#E1, E2 = symbols('E1 E2',real=True, positiv=True)
alpha = symbols('alpha' ,real=True, positiv=True)
S1, S2, kS1, kS2 = symbols('S1 S2 kS1 kS2',real=True, positiv=True)       # potential decomposition fluxes kS
kN1, kN2, km1, km2 = symbols('kN1 kN2 km1 km2',real=True, positiv=True)
aE,B = symbols('aE B',real=True, positiv=True)
cnS1, cnS2, cnE, cnB = symbols('cnS1 cnS2 cnE cnB',real=True, positiv=True)
i2 = symbols('i2',real=True, positiv=True)

Ep = aE*B
Ep = symbols('Ep',real=True, positiv=True)   # rates of enzyme and biomass synthesis and maintenance respiration
E1 = alpha * Ep / kN1
E2 = (1-alpha) * Ep / kN2



#------------ biomass C limited
m, tau, eps = symbols('m tau eps',real=True, positiv=True)   # maintenance and microbial turnover
d1p = kS1*S1
d2p = kS2*S2
d1p, d2p = symbols('d1p d2p',real=True, positiv=True)   # bit simplified
d1 = d1p*E1/(km1+E1)
d2 = d2p*E2/(km2+E2)

#------------ optimal biomass ratio
synE, synB = symbols('synE synB',real=True, positiv=True)   # rates of enzyme and biomass synthesis and maintenance respiration
knB, cnOpt = symbols('knB cnOpt',real=True, positiv=True)   # rates of enzyme and biomass synthesis and maintenance respiration
dE1 = knB*kN1*E1
dE2 = knB*kN2*E2
rMaint = m*B 
#cnOpt = (cnE*synE + cnB*synB)/(synE+synB)
eq1 =  Eq( eps*(d1+d2+dE1+dE2-rMaint)/(d1/cnS1 + d2/cnS2 +dE1/cnE + dE2/cnE) -cnOpt)
#res = solve(eq1, alpha)

eq2 = eq1.subs(kN2, kN1).subs(km2,km1) #.simplify()
#res = solve(eq2, alpha)

eq3 =  Eq( eps*(d1+d2-rMaint)/(d1/cnS1 + d2/cnS2) -cnOpt)
eq3b = eq3.subs(kN2, kN1).subs(km2,km1).simplify()
#res = solve(eq3b, alpha)
