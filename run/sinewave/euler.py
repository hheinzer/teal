from sympy import symbols, sin, Matrix, ccode, simplify

# define independent variables
gamma, x, y, z, t, a, b, c = symbols("gamma x[X] x[Y] x[Z] time a b c")

# define conserved variables (manufactured solution)
rho = 2 + a * sin(b * (x + y + z) - c * t)
rhou = rho
rhov = rho
rhow = rho
rhoe = rho**2

# compute primitive variables
u = rhou / rho
v = rhov / rho
w = rhow / rho
p = (gamma - 1) * (rhoe - rho * (u**2 + v**2 + w**2) / 2)

# compute source terms
W = Matrix([rho, rhou, rhov, rhow, rhoe])
F = Matrix([rhou, rhou * u + p, rhou * v + 0, rhou * w + 0, u * (rhoe + p)])
G = Matrix([rhov, rhov * u + 0, rhov * v + p, rhov * w + 0, v * (rhoe + p)])
H = Matrix([rhow, rhow * u + 0, rhow * v + 0, rhow * w + p, w * (rhoe + p)])
Q = W.diff(t) + F.diff(x) + G.diff(y) + H.diff(z)
for i, var in enumerate(["D", "DU", "DV", "DW", "DE"]):
    print(f"q[{var}] =", ccode(simplify(Q[i])), ";")
