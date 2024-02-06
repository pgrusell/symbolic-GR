'''Calculation of the tensors for the euclidean metrics in
spherical coordinates.'''

import sympy as sy 
from GR_functions import *

sy.init_printing(use_unicode=True)

# Variables
R, theta, phi, x, y = sy.symbols('R, theta, phi, x, y',real=True, positive=True)

# Tranformation 
x = R * sy.cos(phi) * sy.sin(theta)
y = R * sy.sin(phi) * sy.sin(theta)
z = R * sy.cos(theta)

X = sy.Matrix([x, y, z])
Y = sy.Matrix([theta, phi])

# Metrics
g = MetricaTransformada(X, Y, sy.diag(1,1,1))

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