from sympy import *

# define independent variables
gamma, x, y, t, a, b, c = symbols("gamma x y t a b c")

# define manufactured solution (conservative)
rho = 2 + a * sin(b * (x + y) - c * t)
rhou = rho
rhov = rho
rhoe = rho**2

# compute primitive variables
u = rhou / rho
v = rhov / rho
p = (gamma - 1) * (rhoe - rho * (u**2 + v**2) / 2)

# compute source terms
W = Matrix([rho, rho * u, rho * v, rhoe])
F = Matrix([rho * u, rho * u**2 + p, rho * u * v, u * (rhoe + p)])
G = Matrix([rho * v, rho * u * v, rho * v**2 + p, v * (rhoe + p)])
Q = W.diff(t) + F.diff(x) + G.diff(y)
for i in range(4):
    print(f"q[{i}] =", simplify(Q[i]))
