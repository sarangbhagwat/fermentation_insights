# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:48:40 2024

@author: sarangbhagwat
"""
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import contourplots
import itertools
from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap
# from fermentation_insights.plots.analyze_all_combinations import coeff

#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

product_IDs = [
               'TAL', 'TAL_SA',
                'HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol',
               'succinic', 'succinic_neutral']
feedstock_IDs = ['glucose', 'sugarcane', 'corn', 'cornstover']


all_pf_combinations = list(itertools.product(product_IDs, feedstock_IDs))

all_filenames = []
refinery = {i: {} for i in feedstock_IDs}
for p,f in all_pf_combinations:
    all_filenames.append(p+'_'+f)
    refinery[f][p] = np.load(all_filenames[-1]+'_coefficients'+'.npy')
    
# for filename in [f'TAL_SA_{str(i)}x_blp_glucose' for i in [0.1, 0.5, 2, 10]]:
#     all_filenames.append(filename)
#     # refinery[f][p] = np.load(all_filenames[-1])
    
coeff = {i: np.load(i+'_coefficients'+'.npy') for i in all_filenames}


#%% MPSP and titer functions; hyperbola property utils

def MPSP_f(y, t, a, b, c, d):
    return a + b/y + c/t + d/(y*t)

def titer_f(y, MPSP, a, b, c, d):
    inv_MPSP_minus_a = 1/(MPSP-a)
    y_minus_h = y - b*inv_MPSP_minus_a
    y_minus_h[np.where(y_minus_h<=0)]=1e-6
    return A(MPSP, a, b, c, d) / y_minus_h\
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
yields = yields_for_plot = np.linspace(0.01, 200., 100000)
titers = titers_for_plot = np.linspace(0.2, 200, 50)

refinery = 'succinic_glucose'
name = refinery

g = coeff[name][4]
a = coeff[name][0]/g
b = coeff[name][1]/g
c = coeff[name][2]/g
d = coeff[name][3]/g

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

y_label = r"$\bfTiter$" # title of the y axis
y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"

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
#%%
def round_to(val, round_to):
    if isinstance(val, str):
        return val
    else:
        if round_to<1:
            return int(np.round(val,round_to))
        else:
            return np.round(val,round_to)
        
#%%
MPSP = 3

plt.rcParams['font.sans-serif'] = "Arial Unicode"
plt.rcParams['font.size'] = str(12)

full_x_label = x_label + " [" + x_units + "]"
full_y_label = y_label + " [" + y_units + "]" + " at which \nMPSP = " + r"$\mathrm{\$}"+f"{round(MPSP,2)}$"+r"$\cdot\mathrm{kg}^{-1}$"
old_color = 'black'
# mid_color = '#44464a'
mid_color = '#c2c0c0'
new_color = 'blue'
asymptote_linestyle = 'dashed'
asymptote_alpha = 1.
redraw_unchanged_asymptotes = False
panels_together = False

#%%
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
y_ticks = [0, 20, 40, 60, 80, 100, 120]

def format_ax(ax, 
              label_x=True,
              label_y=True, 
              label_top_ticks=True, 
              round_xticks_to=2,
              round_yticks_to=0,
              x_ticks=x_ticks,
              y_ticks=y_ticks,
              axis_title_fontsize=14,
              ):
    if label_x: ax.set_xlabel(full_x_label, fontsize=14)
    if label_y: ax.set_ylabel(full_y_label, fontsize=14)
    ax.tick_params( direction = 'inout' , which='both')
    ax.set_xlim((x_ticks[0], x_ticks[-1]))
    ax.set_ylim((y_ticks[0], y_ticks[-1]))
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    if not label_x: ax.tick_params(labelbottom=False)
    if not label_y: ax.tick_params(labelleft=False)
    
    ax_r = ax.secondary_yaxis('right')
    ax_t = ax.secondary_xaxis('top')
    
    ax_r.tick_params(axis='y', direction='in', which='both')
    ax_r.tick_params(labelright=False)   
    ax_r.set_ylim((y_ticks[0], y_ticks[-1]))
    ax_r.set_yticks(y_ticks)
    
    ax_t.tick_params(axis='x', direction='in', which='both')
    ax_t.tick_params(labeltop=False)   
    ax_t.set_xlim((x_ticks[0], x_ticks[-1]))
    ax_t.set_xticks(x_ticks)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax_r.yaxis.set_minor_locator(AutoMinorLocator())
    ax_t.xaxis.set_minor_locator(AutoMinorLocator())
    if not label_top_ticks:
        # x_ticks = ax.get_xticks()
        # y_ticks = ax.get_yticks()
        xticks_new = [round_to(i,round_xticks_to) for i in x_ticks.copy()]
        xticks_new[-1] = ''
        ax.set_xticklabels(xticks_new)
        yticks_new = [round_to(i,round_yticks_to) for i in y_ticks.copy()]
        yticks_new[-1] = ''
        ax.set_yticklabels(yticks_new)

def mark_asymptotes(ax, h, k, color, zorder=10):
    ax.scatter(h, y_ticks[-1]*1.02, color=[color], marker="v", zorder=zorder, clip_on=False, alpha=asymptote_alpha)
    ax.scatter(x_ticks[-1]*1.02, k, color=[color], marker="<", zorder=zorder, clip_on=False, alpha=asymptote_alpha)

#%% 

fig, axs = plt.subplots(2, 2, constrained_layout=True)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# plt.xlim(0, 1)
# plt.ylim(0, 120)

########### When a changes
ax = axs[0][0]

a_vals = np.linspace(a, 7.5*a, 10)

if panels_together:
    format_ax(ax, label_x=False, label_y=True, label_top_ticks=True,)
else:
    format_ax(ax)

ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a_vals[0], b=b, c=c, d=d), color=old_color)
for a_val in a_vals[1:-1]:
    ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a_val, b=b, c=c, d=d), color=mid_color)
ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a_vals[-1], b=b, c=c, d=d), color=new_color)

h1 = h(MPSP=MPSP, a=a_vals[0], b=b, c=c, d=d)
k1 = k(MPSP=MPSP, a=a_vals[0], b=b, c=c, d=d)
mark_asymptotes(ax, h1, k1, old_color)

h2 = h(MPSP=MPSP, a=a_vals[-1], b=b, c=c, d=d)
k2 = k(MPSP=MPSP, a=a_vals[-1], b=b, c=c, d=d)
mark_asymptotes(ax, h2, k2, new_color, zorder=9)

# plt.show()

######### When b changes
ax = axs[0][1]
b_vals = np.linspace(b, 4*b, 10)

if panels_together:
    format_ax(ax, label_x=False, label_y=True, label_top_ticks=True,)
else:
    format_ax(ax)

ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b_vals[0], c=c, d=d), color=old_color)
for b_val in b_vals[1:-1]:
    ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b_val, c=c, d=d), color=mid_color)
ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b_vals[-1], c=c, d=d), color=new_color)

h1 = h(MPSP=MPSP, a=a, b=b_vals[0], c=c, d=d)
k1 = k(MPSP=MPSP, a=a, b=b_vals[0], c=c, d=d)
mark_asymptotes(ax, h1, k1, old_color)

h2 = h(MPSP=MPSP, a=a, b=b_vals[-1], c=c, d=d)
k2 = k(MPSP=MPSP, a=a, b=b_vals[-1], c=c, d=d)
mark_asymptotes(ax, h2, k2, new_color, zorder=9)

# plt.show()
########### When c changes
ax = axs[1][0]

c_vals = np.linspace(c, 15*c, 10)

if panels_together:
    format_ax(ax, label_x=False, label_y=True, label_top_ticks=True,)
else:
    format_ax(ax)

ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c_vals[0], d=d), color=old_color)
for c_val in c_vals[1:-1]:
    ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c_val, d=d), color=mid_color)
ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c_vals[-1], d=d), color=new_color)

h1 = h(MPSP=MPSP, a=a, b=b, c=c_vals[0], d=d)
k1 = k(MPSP=MPSP, a=a, b=b, c=c_vals[0], d=d)
mark_asymptotes(ax, h1, k1, old_color)

h2 = h(MPSP=MPSP, a=a, b=b, c=c_vals[-1], d=d)
k2 = k(MPSP=MPSP, a=a, b=b, c=c_vals[-1], d=d)
mark_asymptotes(ax, h2, k2, new_color, zorder=9)

# plt.show()
############ When d changes

ax = axs[1][1]

d_vals = np.linspace(d, 30*d, 10)

if panels_together:
    format_ax(ax, label_x=False, label_y=True, label_top_ticks=True,)
else:
    format_ax(ax)

ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c, d=d_vals[0]), color=old_color)
for d_val in d_vals[1:-1]:
    ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c, d=d_val), color=mid_color)
ax.plot(yields, titer_f(yields, MPSP=MPSP, a=a, b=b, c=c, d=d_vals[-1]), color=new_color)

h1 = h(MPSP=MPSP, a=a, b=b, c=c, d=d_vals[0])
k1 = k(MPSP=MPSP, a=a, b=b, c=c, d=d_vals[0])
mark_asymptotes(ax, h1, k1, old_color)

h2 = h(MPSP=MPSP, a=a, b=b, c=c, d=d_vals[-1])
k2 = k(MPSP=MPSP, a=a, b=b, c=c, d=d_vals[-1])
mark_asymptotes(ax, h2, k2, new_color, zorder=9)

# plt.show()
######### final formatting and save
if panels_together: plt.subplots_adjust(wspace=0, hspace=0)
fig.set_figheight(8.5)
fig.set_figwidth(9)

plt.savefig(f'coeffs_description_single_contour.png', 
            transparent = False,  
            facecolor = 'white',
            bbox_inches='tight',
            dpi=600,
            )  

plt.show()

