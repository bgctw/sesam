# derving explicit formulas for optimal enzyme allocation alpha
# Appendix C in Wutzler23_sesam_alloc

from sympy import *
init_printing()
# x, y, z, t = symbols('x y z t')
# k, m, n = symbols('k m n', integer=True)
# f, g, h = symbols('f g h', cls=Function)

alpha_L, alpha_R, alpha_P, d_L, d_R, d_P, e_P = symbols('alpha_L alpha_R alpha_P d_L d_R d_P e_P', nonnegative=True)
k, a, B, d_E = symbols('k a B d_E', positive=True)

du_L = a*B*k*d_L/(k + alpha_L*a*B)**2
du_R = a*B*k*d_R/(k + alpha_R*a*B)**2
du_P = a*B*k*d_P/(k + e_P +alpha_P*a*B)**2

# LR
eq1 = Eq(du_L, du_R.subs(alpha_R, 1-alpha_L))
eq1
s1 = solve(eq1, alpha_L)
s1
print(latex(s1[1]))

# LRP: du_L == du_R
eq2 = Eq(du_L, du_R.subs(alpha_R, 1-alpha_L-alpha_P))
eq2
s2 = solve(eq2, alpha_L)
s2
print(latex(s2[1]))
sAlphaLP = s2[1]
#sAlphaLP = s2[0] # other root

# LRP: du_L(alpa_P) == du_P
eq3 = Eq(du_L.subs(alpha_L, sAlphaLP), du_P)
eq3
s3 = solve(eq3, alpha_P)
s3
print(latex(s3[1]))

# unsuccessful simplifying attempts
# s3_denom = fraction(s3[1])[1]
# simplify(s3_denom)

print(str(s3[1]))

# LP: du_L == du_P
eq4 = Eq(du_L, du_P.subs(alpha_P, 1-alpha_L))
eq4
s4 = solve(eq4, alpha_L)
s4
print(latex(s4[1]))





