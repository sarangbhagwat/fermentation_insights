# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 02:12:23 2024

@author: sarangbhagwat
"""
import numpy as np
from matplotlib import pyplot as plt
from flexsolve import IQ_interpolation
from scipy.optimize import minimize, curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable

ceil = np.ceil
sqrt = np.sqrt
square = np.square
np_sum = np.sum

#%% Geometry utils
def dist(x1, y1, x2, y2):
    return sqrt(square(x1-x2)+square(y1-y2))

def straight_line(x, m, c):
    return m*x + c

#%% Fit utils
def get_Rsq(indicator_orig, indicator_fit):
    y_mean = indicator_orig.mean()
    TSS, RSS = 0, 0
    for yi, y_pred in zip(indicator_orig, indicator_fit):
        RSS += (yi - y_pred)**2
        TSS += (yi - y_mean)**2
    return 1 - RSS/TSS

#%% Linearity utils
def simple_linearity(curve_f, xs):
    ys = curve_f(xs)
    curve_length = np_sum([dist(xs[i], ys[i], xs[i+1], ys[i+1]) for i in range(len(xs)-1)])
    # print(curve_length)
    straight_line_length = dist(xs[0], ys[0], xs[-1], ys[-1])
    # print(straight_line_length)
    return straight_line_length/curve_length

def curve_fit_linearity(curve_f, xs):
    curve_ys = curve_f(xs)
    ((m, c), r) = curve_fit(straight_line, xs, curve_ys, p0=None)
    straight_line_ys = straight_line(xs, m, c)
    return get_Rsq(curve_ys, straight_line_ys)

#%% Cost curve utils
def ic(S, S0, S_max, n, c0=1.):
    N = ceil(S / S_max) if np.all(S>S_max) else 1
    S_each = S/N
    return N * c0 * (S_each/S0)**n

def ic_normalized(S_by_S0, S_max_by_S0, n, c0=1.):
    return ic(S_by_S0, 1., S_max_by_S0, n, c0)

#%%
S_by_S0s = np.linspace(0.001, 10000, 1000)
ns = np.linspace(0.01, 1., 10)

y_bounds = [0.05, 0.95]
t_bounds = [5, 300]
# multiplier_for_linearity = (y_bounds[1]/t_bounds[0])/(y_bounds[0]/t_bounds[1]) # for y/t
multiplier_for_linearity = y_bounds[1] / y_bounds[0] # for y
# additive_term_for_linearity = 500

linearity_f = curve_fit_linearity
# linearity_f = simple_linearity
#!!! TODO: vectorize
ics = []
linearities = []

for n in ns:
    ics.append([])
    linearities.append([])
    for s in S_by_S0s:
        curve_f = lambda s: ic_normalized(S_by_S0=s,
                                          S_max_by_S0=100, 
                                          n=n, 
                                          c0=1000.)
        linearity = linearity_f(curve_f=curve_f,
                                xs=np.linspace(s, 
                                                multiplier_for_linearity*s, 
                                                # s+additive_term_for_linearity,
                                                1000))
        ics[-1].append(curve_f(s))
        linearities[-1].append(linearity)
        
#%% Plot
ax = plt.subplot()
im = ax.contourf(S_by_S0s, ns, 
                 linearities)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
# ax.clabel(im, fmt='%3.0f', colors='black', fontsize=12)

plt.colorbar(im, cax=cax, label=f'Linearity for {np.round(100*(multiplier_for_linearity-1),2)}% increase')

plt.figure()
