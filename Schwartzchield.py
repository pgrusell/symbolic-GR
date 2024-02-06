'''Calculation of the tensors for the Schwartzchield metrics'''

import sympy as sy 
from GR_functions import *

sy.init_printing(use_unicode=True)

# Metrics
r, theta, R, z, phi, t = sy.symbols('r, theta, R, z, phi, t')
variables = sy.Matrix([[t, r, theta, phi]])

A = sy.Function('A')(r)
B = sy.Function('B')(r)
g = sy.diag(-B ,A, r**2, r**2 * sy.sin(theta)**2)

# Calculation
gamma, riemann, ricci, CE = RG(g, variables)

# Display
print('The Christoffel symbols are:')
sy.pprint(gamma)

print('The Riemann tensor is:')
sy.pprint(riemann)

print('The Ricci tensor is:')
sy.pprint(ricci)

print('The curvature scalar is:')
sy.pprint(CE)


