from sympy import symbols, var, simplify, solveset
var('x0,i,k,ts,kc')

Aanalytic =  (x0 - i/k)*1/k*(1-exp(-k*ts)) + i/k*ts
Alinear = ts*(x0 - ts/2*(kc*x0 - i))

s = solveset(Eq(Aanalytic,Alinear), kc)

str(s)
print(latex(s))
