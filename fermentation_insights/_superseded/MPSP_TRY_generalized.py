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
        shifted_rect_hyperbola_f,\
        shifted_rect_hyperbola_two_param,\
        plot_across_param
from biorefineries.TAL._general_utils import add_metrics_to_unit_groups,\
        TEA_breakdown as general_TEA_breakdown
from matplotlib import pyplot as plt
import contourplots

#%%
def create_function(code, namespace_dict):
    def wrapper_fn(statement):
        def f():
            loc={}
            exec(codify(statement), namespace_dict, loc)
            return loc['value_to_return']
        return f
    function = wrapper_fn(code)
    return function

def process_breakdown_across_titers(system, 
                            spec,
                            product,
                            product_ID,
                            TEA_breakdown,
                            titers=np.linspace(20, 90, 20),
                            load_functions=[],
                            simulate_functions=[],
                            recovery_evaluate_function=None,
                            simulate_n_times=3,
                            solve_price_n_times=3,
                            plot=True,
                            exclude_zero_value_metrics_from_plots=True):
    initial_breakdown = TEA_breakdown(False)
    cost_categories = list(initial_breakdown.keys())
    unit_groups = list(initial_breakdown[cost_categories[0]].keys())
    xlabel=f'{product_ID} Titer [g/L]'
    
    evaluate_functions = [lambda: product.cost/product.imass[product_ID]]
    if recovery_evaluate_function: evaluate_functions.append(recovery_evaluate_function)
    cost_categories_begin_index = len(evaluate_functions)
    
    namespace_dict={'TEA_breakdown': TEA_breakdown}
    for ci in cost_categories:
        for ui in unit_groups:
            statement = f"value_to_return = TEA_breakdown(False)['{ci}']['{ui}']"
            evaluate_functions.append(create_function(statement, namespace_dict))
    
    results  = loop_evaluate_across_titers(system=system, spec=spec,
                                product=product, product_ID=product_ID,
                                titers=titers,
                                load_functions=load_functions,
                                evaluate_functions=evaluate_functions,
                                simulate_functions=simulate_functions,
                                simulate_n_times=simulate_n_times,
                                solve_price_n_times=solve_price_n_times,
                                plot=False)
    
    plot_across_param(titers, [results[0]], 
                       ylabel=f'MPSP [$/kg {product_ID}]',
                       xlabel=xlabel)
    
    plot_across_param(titers, [results[1]], 
                       ylabel=f'Recovery of {product_ID} [%]',
                       xlabel=xlabel)
    
    labels = unit_groups
    n_labels = len(labels)
    ylabels = []
    for ci in cost_categories:
        for i in range(n_labels):
            ylabels.append(ci)
    
    labels *= len(cost_categories)
    
    for i in range(len(cost_categories)):
        plot_across_param(titers, results[cost_categories_begin_index+i*n_labels:cost_categories_begin_index+(i+1)*n_labels], 
                           labels=unit_groups, 
                           ylabel=cost_categories[i],
                           xlabel=xlabel,
                           exclude_zero_value_metrics_from_plots=True,)
    
    return results


def process_breakdown_across_yields(system, 
                            spec,
                            product,
                            product_ID,
                            TEA_breakdown,
                            yields=np.linspace(0.2, 0.9, 20),
                            load_functions=[],
                            simulate_functions=[],
                            recovery_evaluate_function=None,
                            simulate_n_times=3,
                            solve_price_n_times=3,
                            plot=True,
                            exclude_zero_value_metrics_from_plots=True):
    initial_breakdown = TEA_breakdown(False)
    cost_categories = list(initial_breakdown.keys())
    unit_groups = list(initial_breakdown[cost_categories[0]].keys())
    xlabel=f'{product_ID} Yield [% theoretical]'
    
    evaluate_functions = [lambda: product.cost/product.imass[product_ID]]
    if recovery_evaluate_function: evaluate_functions.append(recovery_evaluate_function)
    cost_categories_begin_index = len(evaluate_functions)
    
    namespace_dict={'TEA_breakdown': TEA_breakdown}
    for ci in cost_categories:
        for ui in unit_groups:
            statement = f"value_to_return = TEA_breakdown(False)['{ci}']['{ui}']"
            evaluate_functions.append(create_function(statement, namespace_dict))
    
    results  = loop_evaluate_across_yields(system=system, spec=spec,
                                product=product, product_ID=product_ID,
                                yields=yields,
                                load_functions=load_functions,
                                evaluate_functions=evaluate_functions,
                                simulate_functions=simulate_functions,
                                simulate_n_times=simulate_n_times,
                                solve_price_n_times=solve_price_n_times,
                                plot=False)
    
    plot_across_param(yields, [results[0]], 
                       ylabel=f'MPSP [$/kg {product_ID}]',
                       xlabel=xlabel)
    
    plot_across_param(yields, [results[1]], 
                       ylabel=f'Recovery of {product_ID} [%]',
                       xlabel=xlabel)
    
    labels = unit_groups
    n_labels = len(labels)
    ylabels = []
    for ci in cost_categories:
        for i in range(n_labels):
            ylabels.append(ci)
    
    labels *= len(cost_categories)
    
    
    for i in range(len(cost_categories)):
        plot_across_param(yields, results[cost_categories_begin_index+i*n_labels:cost_categories_begin_index+(i+1)*n_labels], 
                           labels=unit_groups, 
                           ylabel=cost_categories[i],
                           xlabel=xlabel,
                           exclude_zero_value_metrics_from_plots=True,)
    
    return results


#%% TAL TRY
from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL import results_metric_1,\
                                                            yields, titers, titers_mol_per_mol_total, productivities,\
                                                            colors, CABBI_green_colormap, get_rounded_str,\
                                                            TAL_maximum_viable_market_range as market_range

#%% HP TRY
# from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
#                                                             yields, titers, titers_mol_per_mol_total, productivities,\
#                                                             colors, CABBI_green_colormap, get_rounded_str,\
#                                                             AA_market_range as market_range,\
#                                                             spec, get_AA_MPSP

#%%
titers_mol_mol = False
titers_gpL = None
if titers_mol_mol:
    titers_gpL = titers.copy()
    titers = np.array(titers_mol_per_mol_total[10])
    
#%% Across titer at alternative yields
product_ID = 'HP'
ms_fit_coeffs, cs_fit_coeffs = {}, {}
indicator_array = np.array(results_metric_1)

#%%
indices = []
ys, ts, mpsps = [], [], []
# find random indices with non-nan results
while len(indices)<4:
    ind = [np.random.randint(0, len(yields)-1), np.random.randint(0, len(titers)-1)]
    if not np.any(np.isnan(indicator_array[0][ind[0], ind[1]])):
        indices.append(ind)
        ys.append(yields[ind[0]])
        ts.append(titers[ind[1]])
        mpsps.append(indicator_array[0][ind[0], ind[1]])

ys, ts, mpsps = np.array(ys), np.array(ts), np.array(mpsps)

#%% 
from scipy.optimize import fsolve

def equations(p):
    a, b, c, d = p
    return (a + b/ys[0] + c/ts[0] + d/(ys[0]*ts[0]) - mpsps[0],
            a + b/ys[1] + c/ts[1] + d/(ys[1]*ts[1]) - mpsps[1],
            a + b/ys[2] + c/ts[2] + d/(ys[2]*ts[2]) - mpsps[2],
            a + b/ys[3] + c/ts[3] + d/(ys[3]*ts[3]) - mpsps[3],
            )

a, b, c, d =  fsolve(equations, (1, 1, 1, 1))

            
#%% Coefficient of determination
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

print(get_Rsq(indicator_array, np.array([[[shifted_rect_hyperbola_two_param(y, t, a, b, c, d) for y in yields] for t in titers]])))

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
    # y_units =r"$\mathrm{mol-product} \cdot \mathrm{mol-product-and-water}^{-1}$"
    # y_ticks = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07]


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
    
    #%% More plot stuff

    fps = 3
    axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
    default_fontsize = 11.
    clabel_fontsize = 9.5
    axis_tick_fontsize = 9.5
    keep_frames = True

    print('\nCreating and saving contour plots ...\n')
    
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
                                    
                                    w_data_vs_x_y_at_multiple_z=[[[shifted_rect_hyperbola_two_param(y,t,a, b, c, d) for y in yields] for t in titers]], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                    
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
    
    m,c,r = fit_straight_line(p_curr, p_curr*arr)
    
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
##

