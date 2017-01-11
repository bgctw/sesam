from __future__ import division 
from sympy import symbols, solve, Eq, diff, sqrt, simplify 

#ER, EL = symbols('ER EL',real=True, positiv=True)
alpha = symbols('alpha' ,real=True, positiv=True)
R, L, kR, kL = symbols('R L kR kL',real=True, positiv=True)       # potential decomposition fluxes kS
kNR, kNL, kmR, kmL = symbols('kNR kNL kmR kmL',real=True, positiv=True)
aE,B = symbols('aE B',real=True, positiv=True)
cnR, cnL, cnE, cnB = symbols('cnR cnL cnE cnB',real=True, positiv=True)
#i2 = symbols('i2',real=True, positiv=True)
imm = symbols('imm',real=True, positiv=True)        # N immobilization 

Ep = aE*B
Ep = symbols('Ep',real=True, positiv=True)   # rates of enzyme synthesis 
ER = alpha * Ep / kNR
EL = (1-alpha) * Ep / kNL



#------------ biomass C limited
m, tau, eps = symbols('m tau eps',real=True, positiv=True)   # maintenance and microbial turnover
decPotR = kR*R
decPotL = kL*L
decPotR, decPotL = symbols('decPotR decPotL',real=True, positiv=True)   # bit simplified
dR = decPotR*ER/(kmR+ER)
dL = decPotL*EL/(kmL+EL)

#------------ optimal enzyme allocation ratio
synE, synB = symbols('synE synB',real=True, positiv=True)   # rates of enzyme and biomass synthesis and maintenance respiration
knB, cnOpt = symbols('knB cnOpt',real=True, positiv=True)   # proportion of enzyme tvr to DOM and biomass CN ratio   
dE1 = knB*kNR*ER
dE2 = knB*kNL*EL
#rMaint = m*B 
#cnOpt = (cnE*synE + cnB*synB)/(synE+synB)
eq1 =  Eq( eps*(dR+dL+dE1+dE2-rMaint)/(dR/cnR + dL/cnL +dE1/cnE + dE2/cnE + imm) -cnOpt)
#res = solve(eq1, alpha)

eq2 = eq1.subs(kNL, kNR).subs(kmL,kmR) #.simplify()
#res = solve(eq2, alpha)

rMaint,E = symbols('rMaint E',real=True, positiv=True)   # proportion of enzyme tvr to DOM and biomass CN ratio
ER  = alpha*E
EL = (1-alpha)*E 
dR = decPotR*ER/(kmR+ER)
dL = decPotL*EL/(kmL+EL)
eq3 =  Eq( eps*(dR+dL-rMaint)/(dR/cnR + dL/cnL +imm) -cnOpt)
#res = solve(eq3, alpha)
eq3b = eq3.subs(kNL, kNR).subs(kmL,kmR).simplify()
#res = solve(eq3b, alpha)
