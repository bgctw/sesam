from __future__ import division 
from sympy import * 
x, y, z, t = symbols('x y z t') 
k, m, n = symbols('k m n', integer=True) 
f, g, h = symbols('f g h', cls=Function) 
E = symbols('E')
alpha = symbols('alpha')
parms_cnB = symbols('parms_cnB')
parms_eps = symbols('parms_eps')
respMaint = symbols('respMaint')
parms_cnS1, parms_cnS2, parms_kS, parms_K = symbols('parms_cnS1 parms_cnS2 parms_kS parms_K')
decC1 = parms_kS * alpha * E / (parms_K + alpha*E)
decC2 = parms_kS * (1-alpha) * E / (parms_K + (1-alpha)*E)
parms_cnB = (parms_eps*(decC1+decC2)-respMaint)/(decC1/parms_cnS1 + decC2/parms_cnS2)
alphaOpt = solve( parms_cnB, alpha)
alphaOpt