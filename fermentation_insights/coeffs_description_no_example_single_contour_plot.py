# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:48:40 2024

@author: sarangbhagwat
"""
import numpy as np
import os
from matplotlib import pyplot as plt
import contourplots
import itertools
from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap
from fermentation_insights.analyze_all_combinations import coeff

#%%

def CABBI_green_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.

    """
    CABBI_colors = (colors.CABBI_orange.RGBn,
                    colors.CABBI_yellow.RGBn,

                    colors.CABBI_green.RGBn,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

#%% MPSP and titer functions; hyperbola property utils

def MPSP_f(y, t, a, b, c, d):
    return a + b/y + c/t + d/(y*t)

def titer_f(y, MPSP, a, b, c, d):
    inv_MPSP_minus_a = 1/(MPSP-a)
    return A(MPSP, a, b, c, d) / (y - b*inv_MPSP_minus_a)\
            + c*inv_MPSP_minus_a

sqrt = np.sqrt
# sqrt_2 = sqrt(2)

def A(MPSP, a, b, c, d):
    inv_MPSP_minus_a = 1/(MPSP-a)
    return (d+b*c*inv_MPSP_minus_a)*inv_MPSP_minus_a

def h(MPSP, a, b, c, d):
    return b/(MPSP-a)

def k(MPSP, a, b, c, d):
    return c/(MPSP-a)

def focus(MPSP, a, b, c, d):
    inv_MPSP_minus_a = 1/(MPSP-a)
    sqrt_2_times_sqrt_A = sqrt(2.*A(MPSP, a, b, c, d))
    return (sqrt_2_times_sqrt_A + b*inv_MPSP_minus_a,
            sqrt_2_times_sqrt_A + c*inv_MPSP_minus_a)


def vertex(MPSP, a, b, c, d):
    inv_MPSP_minus_a = 1/(MPSP-a)
    sqrt_A = sqrt(A(MPSP, a, b, c, d))
    return (sqrt_A + b*inv_MPSP_minus_a,
            sqrt_A + c*inv_MPSP_minus_a)

def linear_eccentricity(MPSP, a, b, c, d):
    return 2*sqrt(A(MPSP, a, b, c, d))


#%%
yields = yields_for_plot = np.linspace(0.01, 5., 100000)
titers = titers_for_plot = np.linspace(0.2, 5, 50)

# refinery = 'TAL_SA_cornstover'
# name = refinery+'_coefficients.npy'

# g = coeff[name][4]
# a = coeff[name][0]/g
# b = coeff[name][1]/g
# c = coeff[name][2]/g
# d = coeff[name][3]/g

a, b, c, d = 1, 1, 1, 1


# indicator_array = [[[MPSP_f(y, t, a, b, c, d) 
#                    for y in yields_for_plot] 
#                    for t in titers_for_plot]]

## %% 
# ########
# # corners of triangular infeasible region
# infeas_ll = 100*yields[0], titers[0]
# infeas_ul = 100*yields[0], titers[-1]

# infeas_ur_titer_index = -1
# top_titer_indicator_array = indicator_array[:, -1, :][0]
# has_infeas_region = True
# try:
#     infeas_ur_yield_index = np.where(np.isnan(top_titer_indicator_array))[0][-1] + 1
#     infeas_ur = 100*yields[infeas_ur_yield_index], titers[infeas_ur_titer_index]
# except:
#     has_infeas_region = False
    

# infeas_region_shape = {
#     # coords as tuple of tuples: (color, zorder),
#     (infeas_ll, infeas_ur, infeas_ul): ('white', 2), # infeasible region smoothing
#     } if has_infeas_region else {}

# ########

#%%
##################
# Plot ticks, labels, other details
# Parameters analyzed across

x_label = r"$\bfYield$" # title of the x axis
# x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
x_units =r"$\mathrm{g} \cdot \mathrm{g}^{-1}$"
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.]

y_label = r"$\bfTiter$" # title of the y axis
y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
y_ticks = [0, 40, 80, 120, 160, 200]

# Metrics
MPSP_w_label = r"$\bfMPSP$" # title of the color axis
MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"

# GWP_w_label = r"$\mathrm{\bfGWP}_{\bf100}$"
GWP_w_label = r"$\mathrm{\bfCarbon}$" + " " + r"$\mathrm{\bfIntensity}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"

FEC_w_label = r"$\bfFEC$" # title of the color axis
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"

AOC_w_label = r"$\bfAOC$" # title of the color axis
AOC_units = r"$\mathrm{MM\$}\cdot\mathrm{y}^{-1}$"

TCI_w_label = r"$\bfTCI$" # title of the color axis
TCI_units = r"$\mathrm{MM\$}$"

Purity_w_label = r"$\bfPurity$" # title of the color axis
Purity_units = "wt " + r"$\mathrm{\%}$"

rel_impact_w_label = r"$\bfRelative$" + " " + r"$\bfImpact$"
rel_impact_units = r"$\mathrm{kg-sugars}\cdot\mathrm{L-broth}^{-1}$"

a_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}^{-1}$"
b_units = r"$\mathrm{\$}\cdot\mathrm{kg-sugars}^{-1}$"
c_units = r"$\mathrm{\$}\cdot\mathrm{L-broth}^{-1}$"
d_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}\cdot\mathrm{L-broth}^{-1}\cdot\mathrm{kg-sugars}^{-1}$"

#%%
fps = 3
axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 11.
clabel_fontsize = 16
axis_tick_fontsize = 16
keep_frames = True

print('\nCreating and saving contour plots ...\n')

get_rounded_str = contourplots.utils.get_rounded_str

#################

MPSP_w_levels = np.arange(0., 8.25, 0.25)
MPSP_cbar_ticks = np.arange(0., 8.1, 1.)
MPSP_w_ticks = [1., 1.25, 1.5, 1.75, 2., 2.5,  3., 4, 5, 6, 8,]
# MPSP_w_ticks = [3, 5, 8]
# MPSP_w_ticks = get_w_ticks(indicator_array, MPSP_w_levels, n_ticks=3)
# MPSP_w_ticks = get_w_ticks_from_percentiles(indicator_array, MPSP_w_levels, 
#                                             percentiles=(25, 50, 75),
#                                             cbar_ticks=MPSP_cbar_ticks)

#%% When a changes

MPSP = 2
plt.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c, d=d), color='black')
focus_1 = focus(MPSP=MPSP, a=a, b=b, c=c, d=d)
vertex_1 = vertex(MPSP=MPSP, a=a, b=b, c=c, d=d)
h1 = h(MPSP=MPSP, a=a, b=b, c=c, d=d)
k1 = k(MPSP=MPSP, a=a, b=b, c=c, d=d)
plt.vlines(h1, 0, titers[-1], colors=['black'], linestyles ='dashed')
plt.hlines(k1, 0, yields[-1], colors=['black'], linestyles ='dashed')
plt.scatter(focus_1[0], focus_1[1])
# plt.scatter(vertex_1[0], vertex_1[1])

plt.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b*2, c=c, d=d), color='blue')
focus_2 = focus(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
vertex_2 = vertex(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
h2 = h(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
k2 = k(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
plt.vlines(h2, 0, titers[-1], colors=['blue'], linestyles ='dashed')
plt.hlines(k2, 0, yields[-1], colors=['blue'], linestyles ='dashed')
plt.scatter(focus_2[0], focus_2[1])
# plt.scatter(vertex_2[0], vertex_2[1])


plt.xlim(0, 5)
plt.ylim(0, 5)
plt.show()

#%% When b changes

MPSP = 2
plt.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c, d=d), color='black')
focus_1 = focus(MPSP=MPSP, a=a, b=b, c=c, d=d)
vertex_1 = vertex(MPSP=MPSP, a=a, b=b, c=c, d=d)
h1 = h(MPSP=MPSP, a=a, b=b, c=c, d=d)
k1 = k(MPSP=MPSP, a=a, b=b, c=c, d=d)
plt.vlines(h1, 0, titers[-1], colors=['black'], linestyles ='dashed')
plt.hlines(k1, 0, yields[-1], colors=['black'], linestyles ='dashed')
# plt.scatter(focus_1[0], focus_1[1])
plt.scatter(vertex_1[0], vertex_1[1])

plt.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b*2, c=c, d=d), color='blue')
focus_2 = focus(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
vertex_2 = vertex(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
h2 = h(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
k2 = k(MPSP=MPSP, a=a, b=b*2, c=c, d=d)
plt.vlines(h2, 0, titers[-1], colors=['blue'], linestyles ='dashed')
plt.hlines(k2, 0, yields[-1], colors=['blue'], linestyles ='dashed')
# plt.scatter(focus_2[0], focus_2[1])
plt.scatter(vertex_2[0], vertex_2[1])


plt.xlim(0, 1)
plt.ylim(0, 200)
plt.show()

