# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 01:30:55 2024

@author: sarangbhagwat
"""
import numpy as np
from matplotlib import pyplot as plt
import sympy as sp

integrate = sp.integrate

#%%

S_by_S0, n, c0 = sp.symbols('S_by_S0 n c0')

ic = c0*(S_by_S0)**n
ic_f = sp.utilities.lambdify((S_by_S0, n, c0 ), ic)

# Partial derivative with respect to S_by_S0
d_ic_d_S_by_S0 = sp.diff(ic, S_by_S0) 

to_integrate = (1+ic**2)**0.5

def linearity_f(curve_f, S_by_S0s, n, c0=1):
    x1, x2 = S_by_S0s[0], S_by_S0s[1]
    curve_length = integrate(to_integrate, (S_by_S0, x2, x1))
    y1, y2 = curve_f(x1, n, c0), curve_f(x2, n, c0)
    straight_line_length = ((x2-x1)**2  + (y2-y1)**2)**0.5
    return curve_length/straight_line_length

#%%
S_by_S0s = np.linspace(0.001, 1000., 10000)
# ics = ic(S_by_S0s, 0.1)