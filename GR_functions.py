'''
This module contains all the functions to perform the GR calculations.
'''

import sympy as sy

## Metric

# Christoffel Symbols
def ChristoffelSimbolo(g,variables):
    '''Returns the Christoffel symbols of the metrics g'''
    n, n = sy.shape(g)
    g_inversa = g**-1
    Gamma_lista = []
    for a in range(n):
        for m in range(n):
            for v in range(n):
                s = 0
                for l in range(n):
                    g_mu = sy.diff(g[l,v], variables[m])
                    g_nu = sy.diff(g[m,l], variables[v])
                    g_lambda=sy.diff(g[m,v],variables[l])
                    G = g_inversa[a, l] * 1/2*(g_mu + g_nu - g_lambda)
                    s += G
                Gamma_lista.append(s)
    Gamma = sy.Array(Gamma_lista, (n,n,n))
    return Gamma

# Riemann tensor
def RiemannTensor(g, variables, Gamma=None, asume_gamma=True):
    '''Returns the Riemann tensor from the Christoffel symbols. 
    If not provided, this funtion calculates them from the metrics g.'''
    if asume_gamma:
        Gamma = ChristoffelSimbolo(g, variables)
    n, n=sy.shape(g)
    g_inversa = g**-1
    Riemann_lista = []
    for m in range(n):
        for v in range(n):
            for a in range(n):
                for b in range(n):
                    k = 0
                    G_a = sy.diff(Gamma[m,v,b], variables[a])
                    G_b = sy.diff(Gamma[m,v,a], variables[b])
                    for s in range(n):
                        k += Gamma[m,s,a]*Gamma[s,v,b] - Gamma[m,s,b]*Gamma[s,v,a]
                    Riemann_lista.append(sy.simplify(G_a - G_b+k))
    Riemann = sy.Array(Riemann_lista, (n,n,n,n))
    return Riemann

# Ricci tensor
def RicciTensor(g, variables, Riemann=None, asume_Riemann=True):
    '''Returns the Ricci tensor from the Riemann symbols. 
    If not provided, this funtion calculates it from the metrics g.'''
    if asume_Riemann:
        Riemann = RiemannTensor(g, variables)
    n, n=sy.shape(g)
    g_inversa = g**-1
    Ricci_lista = []
    for m in range(n):
        for v in range(n):
            k = 0
            for l in range(n):
                k += Riemann[l, m, l, v]
            Ricci_lista.append(k)
    Ricci = sy.Array(Ricci_lista, (n,n))
    return Ricci

# Curvature Scalar
def Escalar(g,variables, Ricci=None, asume_Ricci=True):
    '''Returns the curvature scalar from the Riemann symbols. 
    If not provided, this funtion calculates it from the metrics g.'''
    if asume_Ricci:
        Ricci = RicciTensor(g, variables)
    n, n = sy.shape(g)
    g_inversa = g**-1
    R_escalar = 0
    for l in range(n):
        R_escalar += g_inversa[l, l]*Ricci[l, l]
    return sy.simplify(R_escalar)

# Complete
def RG(g,variables):
    '''Returns the Cristoffel symbols, the Riemann and Ricci tensors,
    and the curvature scalar from the metrics g.'''
    Gamma = ChristoffelSimbolo(g, variables)
    Riemann = RiemannTensor(g, variables, Gamma)
    Ricci = RicciTensor(g, variables, Riemann)
    R = Escalar(g, variables, Ricci)
    return Gamma, Riemann, Ricci, R

## Coordinate tranformations
def MetricaTransformada(X, Y, g0):
    '''This functions performs the transformation of the metrics g0,
    expressed initially on the variables X, in the variables Y.'''
    N = sy.shape(X)[0]
    n = sy.shape(Y)[0]
    Metrica_lista = []
    for m in range(n):
        for v in range(n):
            k = 0
            for i in range(N):
                for j in range(N):
                    k += sy.diff(X[i], Y[m]) * sy.diff(X[j], Y[v]) * g0[i,j]
            Metrica_lista.append(sy.simplify(k))
    Metrica = sy.Array(Metrica_lista, (n,n))
    return Metrica.tomatrix()



