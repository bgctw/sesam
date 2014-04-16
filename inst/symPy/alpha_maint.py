from __future__ import division 
from sympy import symbols, solve, Eq 

E = symbols('E',real=True, positiv=True)        # sum of enzyme levels
alpha = symbols('alpha' ,real=True, positiv=True)
cnOpt, cnS1, cnS2 = symbols('cnOpt cnS1 cnS2',real=True, positiv=True)
parms_eps = symbols('parms_eps',real=True, positiv=True)
respMaint = symbols('respMaint',real=True, positiv=True)
parms_kS, parms_km = symbols('parms_kS parms_km',real=True, positiv=True)
decC1 = parms_kS * alpha * E / (parms_km + alpha*E)
decC2 = parms_kS * (1-alpha) * E / (parms_km + (1-alpha)*E)

decPotS1, decPotS2 = symbols('decPotS1 decPotS2',real=True, positiv=True)
decC1 = decPotS1  * alpha * E / (parms_km + alpha*E)
decC2 = decPotS2  * (1-alpha) * E / (parms_km + (1-alpha)*E)


eq3 = Eq( cnOpt - (parms_eps*(decC1+decC2-respMaint))/(decC1/cnS1 + decC2/cnS2) )
alphaOpt = solve( eq3, alpha)
alphaOpt[0]


# same without maintenance
eq2 = Eq( cnOpt - (parms_eps*(decC1+decC2))/(decC1/cnS1 + decC2/cnS2) )
alphaOpt2 = solve( eq2, alpha)
alphaOpt2[0]

# S2 purely Carbon, no maintenance
eq3 = Eq( cnOpt - (parms_eps*(decC1+decC2))/(decC1/cnS1) )
alphaOpt3 = solve( eq3, alpha)
alphaOpt3[0]

# S2 purely Carbon, with maintenance
eq4 = Eq( cnOpt - (parms_eps*(decC1+decC2-respMaint))/(decC1/cnS1) )
alphaOpt4 = solve( eq4, alpha)
alphaOpt4[0]

# try to include Enzyme turnover
parms_kN, parms_kNB, parms_cnE = symbols('parms_kN parms_kNB parms_cnE',real=True, positiv=True)
eq5 = Eq( cnOpt - (parms_eps*(decC1+decC2+parms_kNB*parms_kN*E-respMaint))/(decC1/cnS1 + decC2/cnS2 + (parms_kNB*parms_kN*E)/parms_cnE) )
alphaOpt5 = solve( eq5, alpha)
alphaOpt5[0]


