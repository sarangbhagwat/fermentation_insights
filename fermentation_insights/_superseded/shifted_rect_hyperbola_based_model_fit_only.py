#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# !!! Package name and short description goes here.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np
from biorefineries.TAL.models.model_utils import codify
from fermentation_insights.utils import loop_evaluate_across_param,\
        loop_evaluate_across_titers, loop_evaluate_across_yields,\
        straight_line_f, fit_straight_line,\
        piecewise_linear_f, fit_piecewise_linear,\
        shifted_rect_hyperbola_f, fit_shifted_rect_hyperbola,\
        shifted_rect_hyperbola_two_param, fit_shifted_rect_hyperbola_two_param,\
        plot_across_param
from biorefineries.TAL._general_utils import add_metrics_to_unit_groups,\
        TEA_breakdown as general_TEA_breakdown
from matplotlib import pyplot as plt
import contourplots
import os

from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap
import chaospy
import itertools

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
ind_ind_for_MPSP_y = 1
ind_ind_for_MPSP_t = 0

run_simulations = True
inflection_product_yields = []

#%% TAL TRY
# if run_simulations:
#     from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_sc import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
#                                                                 yields, titers, titers_mol_per_mol_total, productivities,\
#                                                                 colors, CABBI_green_colormap, get_rounded_str,\
#                                                                 R302, spec,\
#                                                                 TAL_maximum_viable_market_range as market_range
    
    
#     inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
#                                           1-R302.regular_microbe_conversion])
#     np.save('TAL_sugarcane_MPSP', results_metric_1)
#     np.save('TAL_sugarcane_GWP', results_metric_2)
#     np.save('TAL_sugarcane_FEC', results_metric_3)
#     np.save('TAL_sugarcane_AOC', results_metric_5)
#     np.save('TAL_sugarcane_TCI', results_metric_6)
#     np.save('TAL_sugarcane_inflection_product_yields', inflection_product_yields)

        
# else:
#     results_metric_1 = np.load('TAL_sugarcane_MPSP')
#     results_metric_2 = np.load('TAL_sugarcane_GWP')
#     results_metric_3 = np.load('TAL_sugarcane_FEC')
#     results_metric_5 = np.load('TAL_sugarcane_AOC')
#     results_metric_6 = np.load('TAL_sugarcane_TCI')
#     inflection_product_yields = np.load('TAL_sugarcane_inflection_product_yields')
            
            
#%% 3-HP TRY
# if run_simulations:
#     from biorefineries.HP.analyses.fermentation.TRY_analysis_corn_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
#                                                                 yields, titers, productivities,\
#                                                                 colors, CABBI_green_colormap, get_rounded_str,\
#                                                                 AA_market_range as market_range,\
#                                                                 R302, spec, get_AA_MPSP
    
#     inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
#     np.save('HP_corn_MPSP', results_metric_1)
#     np.save('HP_corn_GWP', results_metric_2)
#     np.save('HP_corn_FEC', results_metric_3)
#     np.save('HP_corn_AOC', results_metric_4)
#     np.save('HP_corn_TCI', results_metric_5)
#     np.save('HP_corn_inflection_product_yields', inflection_product_yields)
    
# else:
#     results_metric_1 = np.load('HP_corn_MPSP')
#     results_metric_2 = np.load('HP_corn_GWP')
#     results_metric_3 = np.load('HP_corn_FEC')
#     results_metric_4 = np.load('HP_corn_AOC')
#     results_metric_5 = np.load('HP_corn_TCI')
#     inflection_product_yields = np.load('HP_corn_inflection_product_yields')

#%% 3-HP Neutral TRY
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_cornstover_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_cornstover_MPSP', results_metric_1)
    np.save('HP_neutral_cornstover_GWP', results_metric_2)
    np.save('HP_neutral_cornstover_FEC', results_metric_3)
    np.save('HP_neutral_cornstover_AOC', results_metric_4)
    np.save('HP_neutral_cornstover_TCI', results_metric_5)
    np.save('HP_neutral_cornstover_inflection_product_yields', inflection_product_yields)
    
else:
    results_metric_1 = np.load('HP_neutral_cornstover_MPSP')
    results_metric_2 = np.load('HP_neutral_cornstover_GWP')
    results_metric_3 = np.load('HP_neutral_cornstover_FEC')
    results_metric_4 = np.load('HP_neutral_cornstover_AOC')
    results_metric_5 = np.load('HP_neutral_cornstover_TCI')
    inflection_product_yields = np.load('HP_neutral_cornstover_inflection_product_yields')
    

#%% Succinic TRY
# if run_simulations:
#     from biorefineries.succinic.analyses.TRY_analysis import results_metric_1, results_metric_2, results_metric_3,\
#                                                                 yields, titers, productivities,\
#                                                                 colors, CABBI_green_colormap,\
#                                                                 R302, spec, get_product_MPSP
#                                                                 # SA_price_range as market_range,\
#     market_range = []
#     inflection_product_yields = np.array([1-0.186])
#     np.save('succinic_sugarcane_MPSP', results_metric_1)
#     np.save('succinic_sugarcane_GWP', results_metric_2)
#     np.save('succinic_sugarcane_FEC', results_metric_3)
#     # np.save('succinic_sugarcane_AOC', results_metric_4)
#     # np.save('succinic_sugarcane_TCI', results_metric_5)
#     np.save('succinic_sugarcane_inflection_product_yields', inflection_product_yields)
    
# else:
#     results_metric_1 = np.load('succinic_sugarcane_MPSP')
#     results_metric_2 = np.load('succinic_sugarcane_GWP')
#     results_metric_3 = np.load('succinic_sugarcane_FEC')
#     # results_metric_4 = np.load('succinic_sugarcane_AOC')
#     # results_metric_5 = np.load('succinic_sugarcane_TCI')
#     inflection_product_yields = np.load('succinic_sugarcane_inflection_product_yields')
    

#%% 2,3-BDO TRY
# ## yield-titer indices should be flipped for BDO, as below
# product_ID = 'BDO'
# ind_ind_for_MPSP_y = 0
# ind_ind_for_MPSP_t = 1
# ##

# steps = 40
# yields = np.linspace(0.05, 0.95, steps)
# titers = np.linspace(5., 210., steps)
# productivities = np.array([1.00])

# market_range = []

# inflection_product_yields = [
#     1-0.055-0.02,
#     # 1-0.055,
#     ]

# os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//TRY_results//BDO')

# arr = np.load('BDO_TRY_2021.11.9.18.57'+'.npy')

# # results_metric_1 = arr[:, :, 0, :].transpose()/907.185
# # results_metric_2 = arr[:, :, 1, :].transpose()
# # results_metric_3 = arr[:, :, 2, :].transpose()

# results_metric_1 = arr[:, :, 0, :]/907.185
# results_metric_2 = arr[:, :, 1, :]
# results_metric_3 = arr[:, :, 2, :]

# results_metric_1_new, results_metric_2_new, results_metric_3_new = [], [], []
# for i, j, k in zip(results_metric_1, results_metric_2, results_metric_3):
#     results_metric_1_new.append(i.flatten())
#     results_metric_2_new.append(j.flatten())
#     results_metric_3_new.append(k.flatten())
    
# results_metric_1 = np.array([results_metric_1_new])
# results_metric_2 = np.array([results_metric_2_new])
# results_metric_3 = np.array([results_metric_3_new])

#%%
inflection_product_yields = list(inflection_product_yields)
indices_inflection_product_yields = []
to_add = [i for i in inflection_product_yields]
to_add.sort(reverse=True)


for yi, y in zip (range(len(yields)), yields):
    if not to_add: break
    if y>to_add[-1]: # this point is past the point of inflection
        indices_inflection_product_yields.append(yi-1)
        to_add.pop()

if not yields[-1] in inflection_product_yields: 
    inflection_product_yields.append(yields[-1])
    indices_inflection_product_yields.append(len(yields)-1)
    

#%%
titers_mol_mol = False
titers_gpL = None
if titers_mol_mol:
    titers_gpL = titers.copy()
    titers = np.array(titers_mol_per_mol_total[0])
    
#%%
indicator_array = np.array(results_metric_3)
yield_upper_bound_index_for_eval = indices_inflection_product_yields[0]+1 # only evaluates for yield regimes that occur up to this yield index

indicator_array_for_eval = indicator_array[:, :, :yield_upper_bound_index_for_eval]
yields_for_eval = yields[:yield_upper_bound_index_for_eval]
titers_for_eval = titers

#%%
def MPSP_f(y, t, avals, bvals, cvals, dvals): # use the coefficients solved for the corresponding yield regime
    for ii, i, iy in zip(range(len(indices_inflection_product_yields)),
                         indices_inflection_product_yields, 
                         inflection_product_yields):
        if y<=iy:
            return shifted_rect_hyperbola_two_param(y, t, avals[ii], bvals[ii], cvals[ii], dvals[ii])
    # return shifted_rect_hyperbola_two_param(y, t, avals[-1], bvals[-1], cvals[-1], dvals[-1])

#%%
def get_Rsq(indicator_orig, indicator_fit):
    indicator_orig = indicator_orig.flatten()
    indicator_fit = indicator_fit.flatten()
    y_mean = indicator_orig[np.where(~np.isnan(indicator_orig))].mean()
    TSS, RSS = 0, 0
    for yi, y_pred in zip(indicator_orig, indicator_fit):
        if not np.isnan(yi):
            RSS += (yi - y_pred)**2
            TSS += (yi - y_mean)**2
    return 1 - RSS/TSS


#%% Fit eqn to entire data
# non_nan_indices = np.where(~np.isnan(np.array(indicator_array)[0, :, :yield_upper_bound_index_for_eval].flatten()))
# ytarr = np.array(list(itertools.product(yields[:yield_upper_bound_index_for_eval], titers)))[non_nan_indices, :][0].transpose()
# a_fit, b_fit, c_fit, d_fit, r = fit_shifted_rect_hyperbola_two_param((ytarr[0], ytarr[1]), 
#                                                                            np.array(indicator_array)[:, :, :yield_upper_bound_index_for_eval].flatten()[non_nan_indices],
#                                                                            [1,1,1,1,])

ys_fit, ts_fit, mpsps_fit = [], [], []
for i in range(len(indicator_array_for_eval[0])):
    t_curr = titers[i]
    for j in range(len(indicator_array_for_eval[0][i])):
        y_curr = yields[j]
        mpsp_curr = indicator_array_for_eval[0][i][j]
        if not np.isnan(mpsp_curr):
            ys_fit.append(y_curr)
            ts_fit.append(t_curr)
            mpsps_fit.append(mpsp_curr)

a_fit, b_fit, c_fit, d_fit, r = fit_shifted_rect_hyperbola_two_param((ys_fit, ts_fit), 
                                                                           mpsps_fit,
                                                                           [1,1,1,1,])

print(a_fit, b_fit, c_fit, d_fit)
print(get_Rsq(indicator_array_for_eval, 
              np.array([[[MPSP_f(y, t, [a_fit], [b_fit], [c_fit], [d_fit]) 
                          for y in yields_for_eval] 
                             for t in titers_for_eval]])))


#%% Solve numerically for "g", a', b', c', d', e'
from scipy.optimize import minimize


def solve_for_g(ys, ts, mpsps,
                a, b, c, d,
                print_results=True):
    
    def minimize_obj_f(x):
        g, aprime, bprime, cprime, dprime = x
        return sum([(x[1] +
                      x[2]/ys[i] + 
                      x[3]/ts[i] + 
                      x[4]/(ys[i]*ts[i]))**2
                      for i in range(len(ys))])
    
    cons = (
            # {'type': 'ineq', 'fun': lambda x:  x[0]},

            {'type': 'eq', 'fun': lambda x: x[0]*a - x[1]},
            {'type': 'eq', 'fun': lambda x: x[0]*b - x[2]},
            {'type': 'eq', 'fun': lambda x: x[0]*c - x[3]},
            {'type': 'eq', 'fun': lambda x: x[0]*d - x[4]},
            
            # {'type': 'eq', 'fun': lambda x: x[1] + x[2]/ys[0] + x[3]/ts[0] + x[4]/(ys[0]*ts[0]) + x[5]*ts[0] - x[0]*mpsps[0]},
            # {'type': 'eq', 'fun': lambda x: x[1] + x[2]/ys[1] + x[3]/ts[1] + x[4]/(ys[1]*ts[1]) + x[5]*ts[1] - x[0]*mpsps[1]},
            # {'type': 'eq', 'fun': lambda x: x[1] + x[2]/ys[2] + x[3]/ts[2] + x[4]/(ys[2]*ts[2]) + x[5]*ts[2] - x[0]*mpsps[2]},
            # {'type': 'eq', 'fun': lambda x: x[1] + x[2]/ys[3] + x[3]/ts[3] + x[4]/(ys[3]*ts[3]) + x[5]*ts[3] - x[0]*mpsps[3]},
            
            
            )
    
    res = minimize(minimize_obj_f, 
                   (1, 1, 1, 1, 1), 
                   method='SLSQP', 
                    # bounds=(
                    #     (0, None), # a
                    #     (0, None), # b
                    #     (0, None), # c
                    #     (0, None), # d
                    #     ),
                    constraints=cons,
                    tol=1e-8,
                    options={
                        'maxiter':100,
                        'disp':print_results,
                        'ftol':1e-8,
                        }
                   )
    g, aprime, bprime, cprime, dprime = res.x
    if print_results:
        print(g, aprime, bprime, cprime, dprime)
        print(minimize_obj_f((g, aprime, bprime, cprime, dprime)))
        print(mpsps)
        print([shifted_rect_hyperbola_two_param(yt[0], yt[1], aprime/g, bprime/g, cprime/g, dprime/g) for yt in zip(ys, ts)])
        print('\n')
    try:
        rsq_check = get_Rsq(np.array(mpsps), 
                    np.array([shifted_rect_hyperbola_two_param(yt[0], yt[1], aprime/g, bprime/g, cprime/g, dprime/g) for yt in zip(ys, ts)]))
        if print_results: print(rsq_check)
    except Exception as e:
        print('Could not compute R^2 value: {str(e)}.')
    
    Rsq = get_Rsq(indicator_array_for_eval, np.array([[[MPSP_f(y, t, [aprime/g], [bprime/g], [cprime/g], [dprime/g],) for y in yields_for_eval] for t in titers_for_eval]]))
    if print_results: print(Rsq)
    return Rsq, res.success, g, aprime, bprime, cprime, dprime

Rsq, success, g, aprime, bprime, cprime, dprime = solve_for_g(ys_fit[:2], ts_fit[:2], 
                                                                      # mpsps_fit[:4],
                                                                      [shifted_rect_hyperbola_two_param(yt[0], yt[1], a_fit, b_fit, c_fit, d_fit) 
                                                                       for yt in zip(ys_fit[:2], ts_fit[:2])],
                                                                      a_fit, b_fit, c_fit, d_fit)



#%% Relative impact

import sympy as sp
print('\n\n\n')
y, t = sp.symbols('y t')
a,b,c,d = a_fit, b_fit, c_fit, d_fit

f = a + b/y +c/t + d/(y*t)

# Partial derivative with respect to x
df_dy = sp.diff(f, y) 
print(df_dy)  # Output: 2*x + 2*y

# Partial derivative with respect to y
df_dt = sp.diff(f, t) 
print(df_dt)  # Output: 2*x + 2*y

rel_impact = df_dy/df_dt
print(rel_impact)

rel_impact_f = sp.utilities.lambdify((y,t), rel_impact)

from mpl_toolkits.axes_grid1 import make_axes_locatable
ys_plot = yields[:yield_upper_bound_index_for_eval]
ts_plot = titers
ax = plt.subplot()
im = ax.contourf(ys_plot, ts_plot, np.array([[rel_impact_f(yi, ti) for yi in ys_plot] for ti in ts_plot]))
# ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
# plt.xlabel('Chlorite dismutase initial conc. [mol/L]')
# plt.ylabel('Perchlorate reductase initial conc. [mol/L]')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

# ax.clabel(im, fmt='%3.0f', colors='black', fontsize=12)

plt.colorbar(im, cax=cax, label='Relative impact of yield vs titer')

plt.figure()

#%% Plot stuff
plot_MPSP_y_t = True

if plot_MPSP_y_t:
    # Parameters analyzed across

    x_label = r"$\bfYield$" # title of the x axis
    x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
    x_ticks = [20, 30, 40, 50, 60, 70, 80, 90, 100]

    y_label = r"$\bfTiter$" # title of the y axis
    y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
    y_ticks = [20, 30, 40, 50, 60, 70, 80, 90, 100]


    z_label = r"$\bfProductivity$" # title of the z axis
    z_units =  r"$\mathrm{g} \cdot \mathrm{L}^{-1}  \cdot \mathrm{h}^{-1}$"
    z_ticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

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
    
    fps = 3
    axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
    default_fontsize = 11.
    clabel_fontsize = 9.5
    axis_tick_fontsize = 9.5
    keep_frames = True

    print('\nCreating and saving contour plots ...\n')
    
    get_rounded_str = contourplots.utils.get_rounded_str
    #%% MPSP
    
    # MPSP_w_levels, MPSP_w_ticks, MPSP_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
    MPSP_w_levels = np.arange(0., 4.25, 0.25)
    MPSP_cbar_ticks = np.arange(0., 4.1, 1.)
    MPSP_w_ticks = [1., 1.5, 2., 2.5, 2.75, 3., 3.5, 4]
    # MPSP_w_levels = np.arange(0., 15.5, 0.5)
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_1, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                    x_data=100*yields, # x axis values
                                    # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=MPSP_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=MPSP_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=MPSP_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=MPSP_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='max',
                                    cbar_ticks=MPSP_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='MPSP_y_t_'+product_ID, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    comparison_range=market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 3,
                                    # manual_clabels_regular = {
                                    #     MPSP_w_ticks[5]: (45,28),
                                    #     },
                                    # additional_points ={(73, 62.5):('D', 'w', 6)},
                                    fill_bottom_with_cmap_over_color=False, # for TRY
                                    # bottom_fill_bounds = ((0,0), 
                                    #                       (5,60.),
                                    #                       (95,10.)),
                                    # # zoom_data_scale=5,
                                    # text_boxes = {'>4.00': [(5,5), 'white']},
                                    
                                    # add_shapes = {
                                    #     # coords as tuple of tuples: (color, zorder),
                                    #     ((0,0), (20,200), (1,200)): ('white', 2), # infeasible region smoothing
                                    #     },
                                    units_on_newline = (False, False, False, False), # x,y,z,w
                                    units_opening_brackets = [" (",] * 4,
                                    units_closing_brackets = [")",] * 4,
                                    # label_over_color='white',
                                    )
    
    contourplots.animated_contourplot(
                                    # w_data_vs_x_y_at_multiple_z=[[[MPSP_y_t_f(y,t,ms_fit_coeffs, cs_fit_coeffs) for y in yields] for t in titers]], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                    
                                    w_data_vs_x_y_at_multiple_z=[[[MPSP_f(y, t, avals, bvals, cvals, dvals) for y in yields] for t in titers]], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                    
                                    x_data=100*yields, # x axis values
                                    # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=MPSP_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=MPSP_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=MPSP_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=MPSP_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='max',
                                    cbar_ticks=MPSP_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='MPSP_y_t_fit_'+product_ID, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    comparison_range=market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 3,
                                    # manual_clabels_regular = {
                                    #     MPSP_w_ticks[5]: (45,28),
                                    #     },
                                    # additional_points ={(73, 62.5):('D', 'w', 6)},
                                    fill_bottom_with_cmap_over_color=False, # for TRY
                                    # bottom_fill_bounds = ((0,0), 
                                    #                       (5,60.),
                                    #                       (95,10.)),
                                    # # zoom_data_scale=5,
                                    # text_boxes = {'>4.00': [(5,5), 'white']},
                                    
                                    # add_shapes = {
                                    #     # coords as tuple of tuples: (color, zorder),
                                    #     ((0,0), (20,200), (1,200)): ('white', 2), # infeasible region smoothing
                                    #     },
                                    units_on_newline = (False, False, False, False), # x,y,z,w
                                    units_opening_brackets = [" (",] * 4,
                                    units_closing_brackets = [")",] * 4,
                                    # label_over_color='white',
                                    )

#%% Validate against external data
product_ID = 'product'
print('\nValidate with external data:\n')
print('\n#1: Janković et al., normalized annual cost($/kg) vs titer (wt%)')
## data reported in https://doi.org/10.1002/jctb.7690
# product concentration in fermentation broth (wt%)
p_currs = [[1.001303781,
          2.002607562,
          3.003911343,
          4.005215124,
          4.998696219,
          6],
          [0.492829205,
           1.009126467,
           1.603650587,
           2.002607562]
          ]

# normalized annual cost of separation ($/kg)
arrs = [[0.464788732,
       0.615492958,
       0.666197183,
       0.691549296,
       0.708450704,
       0.71971831],
       [0.335211268,
        0.56,
        0.63943662,
        0.666197183,
        ]
       ]


for p_curr, arr in zip(p_currs, arrs):
    p_curr = np.array(p_curr)
    arr = np.array(arr)
    
    # m,c,r = fit_straight_line(p_curr, p_curr*arr)
    m,c,r = fit_shifted_rect_hyperbola(p_curr, arr)
    
    fig = plt.figure()
    ax = plt.subplot(411)
    
    ax.set_ylabel(f'MPSP [$/kg {product_ID}] * titer[g/L]')
    ax.set_xlabel(f'{product_ID} titer [g/L]')
    ax.scatter(p_curr, [p_curr*arr], label='simulated')
    ax.plot(p_curr, [m*t + c for t in p_curr], label='fit') # MPSP * titer vs titer straight line fit
    
    ax = plt.subplot(412)
    ax.set_ylabel(f'MPSP [$/kg {product_ID}]')
    ax.set_xlabel(f'{product_ID} titer [g/L]')
    ax.scatter(p_curr, arr, label='simulated')
    ax.plot(p_curr, [m + c/t for t in p_curr], label='fit') # resulting MPSP vs titer curve
    plt.legend()
    
    ax = plt.subplot(413)
    ax.set_ylabel(f'Residuals [$/kg {product_ID}]')
    ax.set_xlabel(f'{product_ID} titer [g/L]')
    ax.scatter(p_curr, arr - np.array([m + c/t for t in p_curr]))
    
    ax = plt.subplot(414)
    ax.set_ylabel('Standardized residuals [-]')
    ax.set_xlabel(f'{product_ID} titer [g/L]')
    ax.scatter(p_curr, (arr - np.array([m + c/t for t in p_curr]))/arr)
    plt.show()
    
    print(get_Rsq(np.array(arr), 
                  np.array([shifted_rect_hyperbola_f(p, m, c) for p in p_curr])
                  )
          )

    