from __future__ import division 
from sympy import symbols, solve, Eq 

E = symbols('E',real=True, positiv=True)
alpha = symbols('alpha' ,real=True, positiv=True)
parms_cnB, parms_cnB2 = symbols('parms_cnB parms_cnB2',real=True, positiv=True)
parms_eps = symbols('parms_eps',real=True, positiv=True)
respMaint = symbols('respMaint',real=True, positiv=True)
parms_cnS1, parms_cnS2, parms_kS, parms_K = symbols('parms_cnS1 parms_cnS2 parms_kS parms_K',real=True, positiv=True)
decC1 = parms_kS * alpha * E / (parms_K + alpha*E)
decC2 = parms_kS * (1-alpha) * E / (parms_K + (1-alpha)*E)

decPotS1, decPotS2 = symbols('decPotS1 decPotS2',real=True, positiv=True)
decC1 = decPotS1  * alpha * E / (parms_K + alpha*E)
decC2 = decPotS2  * (1-alpha) * E / (parms_K + (1-alpha)*E)


eq3 = Eq( parms_cnB - (parms_eps*(decC1+decC2-respMaint))/(decC1/parms_cnS1 + decC2/parms_cnS2) )
alphaOpt = solve( eq3, alpha)
alphaOpt[0]

# same without maintenance
eq2 = Eq( parms_cnB - (parms_eps*(decC1+decC2))/(decC1/parms_cnS1 + decC2/parms_cnS2) )
alphaOpt2 = solve( eq2, alpha)
alphaOpt2[0]

# S2 purely Carbon, no maintenance
eq3 = Eq( parms_cnB - (parms_eps*(decC1+decC2))/(decC1/parms_cnS1) )
alphaOpt3 = solve( eq3, alpha)
alphaOpt3[0]

# S2 purely Carbon, with maintenance
eq4 = Eq( parms_cnB - (parms_eps*(decC1+decC2-respMaint))/(decC1/parms_cnS1) )
alphaOpt4 = solve( eq4, alpha)
alphaOpt4[0]


