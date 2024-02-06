'''Calculation of the tensors for the euclidean metrics in
ellipsoidal coordinates.'''

import sympy as sy 
from GR_functions import *

sy.init_printing(use_unicode=True)

# Variables
r, theta, R, z, phi, t = sy.symbols('r, theta, R, z, phi, t')
a = sy.symbols('a')

# Tranformation 
x = (r**2 + a**2)**(1/2) * sy.cos(phi) * sy.sin(theta)
y = (r**2 + a**2)**(1/2) * sy.sin(phi) * sy.sin(theta)
z = r * sy.cos(theta)

X = sy.Matrix([t, x, y, z])
Y = sy.Matrix([t, r, theta, phi])

# Metrics
g = MetricaTransformada(X, Y, sy.diag(-1,1,1,1))

# Tensors
gamma, riemann, ricci, CE = RG(g, Y)

# Display
print('The transformed metrics is:')
sy.pprint(g)

print('The Christoffel symbols are:')
sy.pprint(gamma)

print('The Riemann tensor is:')
sy.pprint(riemann)

print('The Ricci tensor is:')
sy.pprint(ricci)

print('The curvature scalar is:')
sy.pprint(CE)