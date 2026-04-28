# manufactured solution for a trigonometric Navier-Stokes flow
from sympy import Matrix, ccode, simplify, sin, symbols

gamma, mu, prandtl, x, y, z, t, a, b, c = symbols("gamma mu prandtl x[X] x[Y] x[Z] time a b c")

# conserved fields
rho = 2 + a * sin(b * (x + y + z) - c * t)
rhou = rho
rhov = rho
rhow = rho
rhoe = rho**2

# primitive fields
u = rhou / rho
v = rhov / rho
w = rhow / rho
p = (gamma - 1) * (rhoe - rho * (u**2 + v**2 + w**2) / 2)

# stress tensor
divU = u.diff(x) + v.diff(y) + w.diff(z)
tauXX = 2 * mu * (u.diff(x) - divU / 3)
tauYY = 2 * mu * (v.diff(y) - divU / 3)
tauZZ = 2 * mu * (w.diff(z) - divU / 3)
tauXY = mu * (u.diff(y) + v.diff(x))
tauXZ = mu * (u.diff(z) + w.diff(x))
tauYZ = mu * (v.diff(z) + w.diff(y))

# heat flux
factor = -mu * gamma / ((gamma - 1) * prandtl * rho**2)
qX = factor * (p.diff(x) * rho - p * rho.diff(x))
qY = factor * (p.diff(y) * rho - p * rho.diff(y))
qZ = factor * (p.diff(z) * rho - p * rho.diff(z))

# state and flux vectors
W = Matrix([rho, rhou, rhov, rhow, rhoe])
FC = Matrix([rhou, rhou * u + p, rhou * v + 0, rhou * w + 0, u * (rhoe + p)])
GC = Matrix([rhov, rhov * u + 0, rhov * v + p, rhov * w + 0, v * (rhoe + p)])
HC = Matrix([rhow, rhow * u + 0, rhow * v + 0, rhow * w + p, w * (rhoe + p)])
FV = Matrix([0, tauXX, tauXY, tauXZ, u * tauXX + v * tauXY + w * tauXZ - qX])
GV = Matrix([0, tauXY, tauYY, tauYZ, u * tauXY + v * tauYY + w * tauYZ - qY])
HV = Matrix([0, tauXZ, tauYZ, tauZZ, u * tauXZ + v * tauYZ + w * tauZZ - qZ])

# source term for the manufactured solution
Q = W.diff(t) + FC.diff(x) + GC.diff(y) + HC.diff(z) - FV.diff(x) - GV.diff(y) - HV.diff(z)
for i, var in enumerate(["D", "DU", "DV", "DW", "DE"]):
    print(f"q[{var}] =", ccode(simplify(Q[i])), ";")
