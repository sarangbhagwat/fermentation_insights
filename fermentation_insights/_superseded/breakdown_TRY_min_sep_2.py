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
        straight_line_f, fit_min_sep_energy_based_f,\
        piecewise_linear_f, fit_piecewise_linear,\
        min_sep_energy_based_f,\
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
# from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL import results_metric_1, yields, titers, productivities,\
#                                                             colors, CABBI_green_colormap, get_rounded_str,\
#                                                             TAL_maximum_viable_market_range

#%% HP TRY
from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                            yields, titers, titers_mol_per_mol_total, productivities,\
                                                            colors, CABBI_green_colormap, get_rounded_str,\
                                                            AA_market_range,\
                                                            spec, get_AA_MPSP

#%% Across yield at alternative titers
# product_ID = 'TAL'
# plot_MPSPs = True
# print_coeffs = True
# plot_coefficients_vs_titer = True
# plot_MPSP_y_t = True
# ms, cs = [], []
# ms_fit_coeffs, cs_fit_coeffs = {}, {}

# for i in range(len(results_metric_1)):
#     p = productivities[i]

#     for j in range(len(results_metric_1[0][0])):
#         t = titers[j]
        
#         where_nan = np.where(np.isnan(results_metric_1[0][j]))
#         from_index = np.max(where_nan) + 1 if list(where_nan[0]) else 0 # 1 + highest yield index at which MPSP=nan 
#                                                                         # (i.e., yield is high enough from this index onwards)
#         MPSPs = results_metric_1[0][j][from_index:]
#         yields_curr = yields[from_index:]
        
#         m,c,r = fit_straight_line(yields_curr, yields_curr*MPSPs)
#         ms.append(m)
#         cs.append(c)
#         if print_coeffs: print(t, ': ', m, c)
        
#         if plot_MPSPs:
#             fig = plt.figure()
#             ax = plt.subplot(211)
            
#             ax.set_ylabel(f'MPSP [$/kg {product_ID}] * yield[% theoretical]')
#             ax.set_xlabel(f'{product_ID} yield [% theoretical]')
#             ax.scatter(yields_curr, [yields_curr*MPSPs], label='simulated')
#             ax.plot(yields_curr, [m*t + c for t in yields_curr], label='fit') # MPSP * yield vs yield straight line fit
            
#             ax = plt.subplot(212)
#             ax.set_ylabel(f'MPSP [$/kg {product_ID}]')
#             ax.set_xlabel(f'{product_ID} yield [% theoretical]')
#             ax.scatter(yields_curr, MPSPs, label='simulated')
#             ax.plot(yields_curr, [m + c/t for t in yields_curr], label='fit') # resulting MPSP vs yield curve
#             plt.legend()
#             plt.show()
            
#     ms, cs = np.array(ms), np.array(cs)
    
#     # Piecewise linear fit for m vs titer
#     ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2'], ms_fit_coeffs['r'] = fit_piecewise_linear(titers, ms*titers)
#     cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2'], cs_fit_coeffs['r'] = fit_piecewise_linear(titers, cs*titers)
    
#     if plot_coefficients_vs_titer:
#         fig = plt.figure()
#         ax = plt.subplot(211)
        
#         ax.set_ylabel(f'm * titer[g/L]')
#         ax.set_xlabel(f'{product_ID} titer[g/L]')
#         ax.scatter(titers, ms*titers, label='orig') # from MPSP*yield vs yield straight line fit
#         ax.plot(titers, piecewise_linear_f(titers, ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2']),
#                 label='fit') # from (orig)*titer vs titer piecewise linear fit
        
#         ax = plt.subplot(212)
#         ax.set_ylabel(f'c * titer[g/L]')
#         ax.set_xlabel(f'{product_ID} titer[g/L]')
#         ax.scatter(titers, cs*titers, label='orig') # from MPSP*yield vs yield straight line fit
#         ax.plot(titers, piecewise_linear_f(titers, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2']),
#                 label='fit') # from (orig)*titer vs titer piecewise linear fit
        
#         plt.legend()
#         plt.show()


# def MPSP_y_t_f(y, t, ms_fit_coeffs, cs_fit_coeffs):
#     t_inv = 1./t
#     m = piecewise_linear_f(t, ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2']) * t_inv
#     c = piecewise_linear_f(t, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2']) * t_inv
#     return shifted_rect_hyperbola_f(y, m, c)

#%%
titers_gpL = titers.copy()
titers = np.array(titers_mol_per_mol_total[0])

#%% Across titer at alternative yields
product_ID = 'HP'
plot_MPSPs = True
print_coeffs = True
plot_coefficients_vs_yield = True
plot_MPSP_y_t = True
# ms, cs = [], []
ms_fit_coeffs, cs_fit_coeffs = {}, {}
indicator_array = np.array(results_metric_1)
# indicator_array = np.array(results_metric_5/(yields*results_metric_6))
fit_or_solve = 'fit'
indices_for_solution = [3,4]

for i in range(len(indicator_array)):
    p = productivities[i]
    ms, cs = [], []
    for j in range(len(indicator_array[0,0,:])):
        y = yields[j]
        
        where_nan = np.where(np.isnan(indicator_array[0,:,j]))
        to_index = np.min(where_nan) if list(where_nan[0]) else -1 # lowest titer index at which MPSP=nan 
                                                                        # (i.e., titer is low enough at lower indices than this)
        MPSPs = indicator_array[0,:,j][:to_index]
        titers_curr = titers[:to_index]
        
        m,c,r = None,None,None
        
        if fit_or_solve=='fit':
            m,c,r = fit_min_sep_energy_based_f(titers_curr, MPSPs)
        elif fit_or_solve=='solve':
            ts = titers_curr[indices_for_solution]
            line_ys = MPSPs[indices_for_solution]
            # m = line_ys[0] - c/ts[0] 
            # line_ys[0] - c/ts[0] + c/ts[1] = line_ys[1]
            c = (ts[0]*ts[1] * (line_ys[1] - line_ys[0]))/(ts[0]-ts[1])
            m = line_ys[0] - c/ts[0] 
            
        ms.append(m)
        cs.append(c)
        if print_coeffs: print(y, ': ', m, c)
        
        if plot_MPSPs:
            fig = plt.figure()
            ax = plt.subplot(411)
            
            ax.set_ylabel(f'MPSP [$/kg {product_ID}] * titer[g/L]')
            ax.set_xlabel(f'{product_ID} titer [g/L]')
            ax.scatter(titers_curr, [titers_curr*MPSPs], label='simulated')
            ax.plot(titers_curr, [m*t + c for t in titers_curr], label='fit') # MPSP * titer vs titer straight line fit
            
            ax = plt.subplot(412)
            ax.set_ylabel(f'MPSP [$/kg {product_ID}]')
            ax.set_xlabel(f'{product_ID} titer [g/L]')
            ax.scatter(titers_curr, MPSPs, label='simulated')
            ax.plot(titers_curr, [min_sep_energy_based_f(t,m,c) for t in titers_curr], label='fit') # resulting MPSP vs titer curve
            plt.legend()
            
            ax = plt.subplot(413)
            ax.set_ylabel(f'Residuals [$/kg {product_ID}]')
            ax.set_xlabel(f'{product_ID} titer [g/L]')
            ax.scatter(titers_curr, MPSPs - np.array([min_sep_energy_based_f(t,m,c) for t in titers_curr]))
            
            ax = plt.subplot(414)
            ax.set_ylabel('Standardized residuals [-]')
            ax.set_xlabel(f'{product_ID} titer [g/L]')
            ax.scatter(titers_curr, (MPSPs - np.array([min_sep_energy_based_f(t,m,c) for t in titers_curr]))/MPSPs)
            plt.show()
            
            
    ms, cs = np.array(ms), np.array(cs)
    
    # Piecewise linear fit for m vs yield
    ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2'], ms_fit_coeffs['r'] = fit_piecewise_linear(yields, ms*yields, p0=[0.5,0.5,1.,1.])
    cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2'], cs_fit_coeffs['r'] = fit_piecewise_linear(yields, cs*yields, p0=[0.5,0.5,1.,1.])
    
    if plot_coefficients_vs_yield:
        fig = plt.figure()
        
        # coeff*yield vs yield
        ax = plt.subplot(211)
        
        ax.set_ylabel(f'm * yield[% theoretical]')
        ax.set_xlabel(f'{product_ID} yield[% theoretical]')
        ax.scatter(yields, ms*yields, label='orig') # from MPSP*titer vs titer straight line fit
        ax.plot(yields, piecewise_linear_f(yields, ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2']),
                label='fit') # from (orig)*yield vs yield piecewise linear fit
        
        ax = plt.subplot(212)
        ax.set_ylabel(f'c * yield[% theoretical]')
        ax.set_xlabel(f'{product_ID} yield[% theoretical]')
        ax.scatter(yields, cs*yields, label='orig') # from MPSP*titer vs titer straight line fit
        ax.plot(yields, piecewise_linear_f(yields, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2']),
                label='fit') # from (orig)*yield vs yield piecewise linear fit
        
        plt.legend()
        plt.show()

def MPSP_y_t_f(y, t, ms_fit_coeffs, cs_fit_coeffs):
    y_inv = 1./y
    yield_index = np.where(yields==y)[0][0]
    titer_index = np.where(titers==t)[0][0]
    m = piecewise_linear_f(y, ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2']) * y_inv
    # m = ms[yield_index]
    c = piecewise_linear_f(y, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2']) * y_inv
    # c = cs[yield_index]
    if np.isnan(results_metric_1[0, titer_index, yield_index]): return np.nan
    return min_sep_energy_based_f(t, m, c)

def mechanistic_MPSP_titer_f(t, a1, a2, t1, n, b):
    return a1 + a2*(t1/t)**n + b/t


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

print(get_Rsq(indicator_array, np.array([[[MPSP_y_t_f(y,t,ms_fit_coeffs, cs_fit_coeffs) for y in yields] for t in titers]])))

#%% Plot stuff

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
                                    comparison_range=AA_market_range,
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
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=[[[MPSP_y_t_f(y,t,ms_fit_coeffs, cs_fit_coeffs) for y in yields] for t in titers]], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
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
                                    comparison_range=AA_market_range,
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



#%%
stop
titers = np.linspace(20, 90, 5)
yields = np.linspace(0.2, 0.5, 5)

#%% TAL

from biorefineries import TAL

TAL.load_TAL_model('A')
TAL.simulate_and_print()
system = TAL.system
spec = TAL.spec

spec.load_specifications(yields[0], spec.baseline_titer)
spec.decrease_byproduct_yields_uniformly = True
# spec.baseline_titer = 18.

product = TAL.system.flowsheet('TAL_product')
product_ID = 'TAL'
TAL_TEA_breakdown = TAL.TEA_breakdown

broth = spec.reactor.outs[1]
results_TAL_titer = process_breakdown_across_titers(system, spec, product, product_ID, TAL_TEA_breakdown, titers,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol[product_ID],
                                                    )


spec.load_specifications(yields[0], spec.baseline_titer, spec.baseline_productivity)
for i in range(3):
    system.simulate()
for i in range(3):
    TAL.get_adjusted_MSP()

results_TAL_yield = process_breakdown_across_yields(system, spec, product, product_ID, TAL_TEA_breakdown, yields,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol[product_ID],
                                                    )


spec.decrease_byproduct_yields_uniformly = False

results_TAL_titer_nyc = process_breakdown_across_titers(system, spec, product, product_ID, TAL_TEA_breakdown, titers,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol[product_ID],
                                                    )

spec.load_specifications(yields[0], spec.baseline_titer, spec.baseline_productivity)
for i in range(3):
    system.simulate()
for i in range(3):
    TAL.get_adjusted_MSP()
    
results_TAL_yield_nyc = process_breakdown_across_yields(system, spec, product, product_ID, TAL_TEA_breakdown, yields,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol[product_ID],
                                                    )

#%% HP

from biorefineries.HP.system_light_lle_vacuum_distillation import HP_sys, HP_tea, HP_lca,\
        spec, AA, process_groups, process_groups_dict

system = HP_sys
spec = spec
product = AA
product_ID = 'AA'

broth = spec.reactor.outs[0]

add_metrics_to_unit_groups(unit_groups=process_groups, system=HP_sys, TEA=HP_tea, LCA=HP_lca)

def HP_TEA_breakdown(print_output=False):
    return general_TEA_breakdown(process_groups_dict, False)

results_HP_titer = process_breakdown_across_titers(system, spec, product, product_ID, HP_TEA_breakdown, titers,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol['HP'],
                                                    )
results_HP_yield = process_breakdown_across_yields(system, spec, product, product_ID, HP_TEA_breakdown, yields,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol['HP'],
                                                    )

#%% Succinic

from biorefineries.succinic.system_sc import succinic_sys, succinic_tea, succinic_LCA,\
        spec, product_stream, unit_groups_dict, TEA_breakdown as succinic_TEA_breakdown

system = succinic_sys
spec = spec
product = product_stream
product_ID = 'SuccinicAcid'

results_succinic_titer = process_breakdown_across_titers(system, spec, product, product_ID, succinic_TEA_breakdown, titers,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol[product_ID],
                                                    )
results_succinic_yield = process_breakdown_across_yields(system, spec, product, product_ID, succinic_TEA_breakdown, yields,
                                                    recovery_evaluate_function=lambda: product.imol[product_ID]/broth.imol[product_ID],
                                                    )

#%% Titer

MPSPs_all = [results_TAL_titer[0], 
             results_TAL_titer_nyc[0],
             # results_HP_titer[0], 
             # results_succinic_titer[0]
             ]
product_IDs = ['TAL', 'TAL', 'AA', 'SuccinicAcid']

# All - MPSP * titer straight line fit
for MPSPs, product_ID in zip(MPSPs_all, product_IDs):
    m,c,r = fit_straight_line(titers, titers*MPSPs)
    fig = plt.figure()
    ax = plt.subplot(111)
    
    ax.set_ylabel(f'MPSP [$/kg {product_ID}] * titer[g/L]')
    ax.set_xlabel(f'{product_ID} titer [% theoretical]')
    ax.scatter(titers, [titers*MPSPs], label='simulated')
    ax.plot(titers, [m*t + c for t in titers], label='fit')
    plt.legend()
    plt.show()
    
# All - MPSP from straight line fit
for MPSPs, product_ID in zip(MPSPs_all, product_IDs):
    m,c,r = fit_straight_line(titers, titers*MPSPs)
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_ylabel(f'MPSP [$/kg {product_ID}]')
    ax.set_xlabel(f'{product_ID} titer [% theoretical]')
    ax.scatter(titers, MPSPs, label='simulated')
    ax.plot(titers, [m + c/t for t in titers], label='fit')
    plt.legend()
    plt.show()

#%% Yield

# All - MPSP * yield straight line fit
MPSPs_all = [results_TAL_yield[0][2:], 
             results_TAL_yield_nyc[0][2:],
             # results_HP_yield[0], 
             # results_succinic_yield[0]
             ]
product_IDs = ['TAL', 'TAL', 'AA', 'SuccinicAcid']
for MPSPs, product_ID in zip(MPSPs_all, product_IDs):
    m,c,r = fit_straight_line(yields[2:], yields[2:]*MPSPs)
    fig = plt.figure()
    ax = plt.subplot(111)
    
    ax.set_ylabel(f'MPSP [$/kg {product_ID}] * yield[g/L]')
    ax.set_xlabel(f'{product_ID} Yield [% theoretical]')
    ax.scatter(yields[2:], [yields[2:]*MPSPs], label='simulated')
    ax.plot(yields[2:], [m*t + c for t in yields[2:]], label='fit')
    plt.legend()
    plt.show()
    
# All - MPSP from straight line fit
for MPSPs, product_ID in zip(MPSPs_all, product_IDs):
    m,c,r = fit_straight_line(yields[2:], yields[2:]*MPSPs)
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_ylabel(f'MPSP [$/kg {product_ID}]')
    ax.set_xlabel(f'{product_ID} Yield [% theoretical]')
    ax.scatter(yields[2:], MPSPs, label='simulated')
    ax.plot(yields[2:], [m + c/t for t in yields[2:]], label='fit')
    plt.legend()
    plt.show()
    

#%% All - MPSP from straight line fit + recovery adjustment
# MPSPs_all = [results_TAL[0], results_HP[0], results_succinic[0]]
# product_IDs = ['TAL', 'AA', 'SuccinicAcid']
# for MPSPs, product_ID in zip(MPSPs_all, product_IDs):
#     m,c,r = fit_straight_line(yields, yields*MPSPs)
#     plot_across_param(yields, [MPSPs, [m + c/t for t in yields]], 
#                        ylabel=f'MPSP [$/kg {product_ID}]',labels=['simulated','straight line fit'])  
