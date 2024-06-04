from numpy import pi, sin, cos, sqrt

# https://homepages.dias.ie/jmackey/jmac/node10.html
# https://en.wikipedia.org/wiki/Moving_shock
# https://en.wikipedia.org/wiki/Thermodynamic_relations_across_normal_shocks
gamma = 1.4
rho0 = 1.4
u0 = 0.0
p0 = 1.0
Mx = 10.0
alpha = 30.0 * pi / 180.0

# compute post-shock speed of sound
a0 = sqrt(gamma * p0 / rho0)

# compute pre-shock state
rho1 = rho0 / (1 - 2 / (gamma + 1) * (1 - 1 / Mx**2))
p1 = p0 * (1 + 2 * gamma / (gamma + 1) * (Mx**2 - 1))
a1 = sqrt(gamma * p1 / rho1)

# compute pre-shock velocity
W = Mx * a0
My = sqrt(((gamma - 1) * Mx**2 + 2) / (2 * gamma * Mx**2 - (gamma - 1)))
uy = My * a1
u1 = W + u0 - uy

# compute rotated pre_shock state
u1x = u1 * cos(alpha)
u1y = u1 * -sin(alpha)
Wx = W * cos(alpha)
Wy = W * -sin(alpha)
print(list(map(lambda x: round(x, 7), [rho1, u1x, u1y, p1, Wx, Wy])))
