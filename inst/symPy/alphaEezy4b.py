from __future__ import division 
from sympy import symbols, solve 

#E1, E2 = symbols('E1 E2',real=True, positiv=True)
alpha = symbols('alpha' ,real=True, positiv=True)
d1, d2 = symbols('d1 d2',real=True, positiv=True)
tau1, tau2, m1, m2 = symbols('tau1 tau2 m1 m2',real=True, positiv=True)
a, B = symbols('a B',real=True, positiv=True)
cn1, cn2, cnE = symbols('cn1 cn2 cnE',real=True, positiv=True)

E1 = alpha * a * B / tau1
E2 = (1-alpha) * a * B / tau2

# C-limited
e1 = d1/(tau1*(E1+m1))      # efficiencies for C acquiring: decomp/E_tvr, E1 cancels
e2 = d2/(tau2*(E2+m2))
# 

# N-limited
e1 = d1/(tau1*(E1+m1)) *  cnE/cn1   
e2 = d2/(tau2*(E2+m2)) *  cnE/cn2    
# (2*B*a*cn2*d1 + cn1*d2*m1*tau1 + cn2*d1*m2*tau2 - sqrt(4*B**2*a**2*cn1*cn2*d1*d2 + 4*B*a*cn1*cn2*d1*d2*m1*tau1 + 4*B*a*cn1*cn2*d1*d2*m2*tau2 + cn1**2*d2**2*m1**2*tau1**2 + 2*cn1*cn2*d1*d2*m1*m2*tau1*tau2 + cn2**2*d1**2*m2**2*tau2**2))/(-2*B*a*cn1*d2 + 2*B*a*cn2*d1)


eq1 = Eq( alpha - e1/(e1+e2))
res = solve(eq1, alpha)
res[0]






