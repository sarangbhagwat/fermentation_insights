#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np
from math import log as ln
import scipy
from matplotlib import pyplot as plt

#%%
def loop_evaluate_across_param(
                          load_functions,
                          evaluate_functions,
                          param_values_to_load,
                          simulate_functions=[lambda: None], # these are called between loading and evaluating
                          plot=True,
                          labels=None,
                          ylabel=None,
                          ):
    results = [[] for i in range(len(evaluate_functions))]
        
    for i in param_values_to_load:
        for lf in load_functions: lf(i)
        for sf in simulate_functions: sf()
        for (ei,ef) in zip(range(len(evaluate_functions)), evaluate_functions): 
            results[ei].append(ef())
            
    results = np.array(results)
    
    if plot: plot_across_param(param_values_to_load, results, labels=labels, ylabel=ylabel)
    
    return results


def loop_evaluate_across_titers(system, 
                                spec,
                                product,
                                product_ID,
                                titers=np.linspace(20, 90, 20),
                                load_functions=[],
                                evaluate_functions=[],
                                simulate_functions=[],
                                simulate_n_times=3,
                                solve_price_n_times=3,
                                plot=True,
                                labels=None,
                                ylabel=None,
                                ):
    TEA = system.TEA
    
    def simulate_system():
        for i in range(simulate_n_times): system.simulate()
    def solve_price():
        for i in range(solve_price_n_times): product.price=TEA.solve_price(product)
        
    if not simulate_functions:
        simulate_functions = [simulate_system, solve_price]
        
    if not evaluate_functions:
        evaluate_functions = [lambda: product.cost/product.imass[product_ID],
                              lambda: TEA.AOC,
                              lambda: TEA.TCI,
                              lambda: TEA.material_cost,
                              lambda: TEA.utility_cost,
                              lambda: product.imass[product_ID]]
        
    if not load_functions:
        load_functions = [lambda t: spec.load_specifications(spec.baseline_yield, t, spec.baseline_productivity)]
    
    return loop_evaluate_across_param(load_functions=load_functions, 
                                      evaluate_functions=evaluate_functions, 
                                      param_values_to_load=titers,
                                      simulate_functions=simulate_functions,
                                      plot=plot,
                                      labels=labels,
                                      ylabel=ylabel)

def loop_evaluate_across_yields(system, 
                                spec,
                                product,
                                product_ID,
                                yields=np.linspace(20, 90, 20),
                                load_functions=[],
                                evaluate_functions=[],
                                simulate_functions=[],
                                simulate_n_times=3,
                                solve_price_n_times=3,
                                plot=True,
                                labels=None,
                                ylabel=None,
                                ):
    TEA = system.TEA
    
    def simulate_system():
        for i in range(simulate_n_times): system.simulate()
    def solve_price():
        for i in range(solve_price_n_times): product.price=TEA.solve_price(product)
        
    if not simulate_functions:
        simulate_functions = [simulate_system, solve_price]
        
    if not evaluate_functions:
        evaluate_functions = [lambda: product.cost/product.imass[product_ID],
                              lambda: TEA.AOC,
                              lambda: TEA.TCI,
                              lambda: TEA.material_cost,
                              lambda: TEA.utility_cost,
                              lambda: product.imass[product_ID]]
        
    if not load_functions:
        load_functions = [lambda y: spec.load_specifications(y, spec.baseline_titer, spec.baseline_productivity)]
    
    return loop_evaluate_across_param(load_functions=load_functions, 
                                      evaluate_functions=evaluate_functions, 
                                      param_values_to_load=yields,
                                      simulate_functions=simulate_functions,
                                      plot=plot,
                                      labels=labels,
                                      ylabel=ylabel)


def plot_across_param(param_values, results, labels=None, ylabel=None, xlabel=None,
                      colors=None, markers=None,
                       exclude_zero_value_metrics_from_plots=False):
    fig = plt.figure()
    ax = plt.subplot(111)
    n_metrics = len(results)
    if not labels: labels = ['' for i in range(n_metrics)]
    for i in range(n_metrics):
        results_to_plot = results[i]
        if not (exclude_zero_value_metrics_from_plots and np.all(results_to_plot==0.)):
            ax.plot(param_values, results_to_plot, label=labels[i], colors=colors, markers=markers)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
#%%

def straight_line_f(x, m, c):
    return m*x + c

def fit_straight_line(xdata, ydata, p0=None):
    ((m,c), r) = scipy.optimize.curve_fit(straight_line_f, xdata, ydata, p0=p0)
    return m, c, r


def piecewise_linear_f(x, x0, y0, m1, m2):
    y = np.piecewise(x, [x < x0, x >= x0],
                     [lambda x:m1*x + y0-m1*x0, lambda x:m2*x + y0-m2*x0])
    return y

def fit_piecewise_linear(xdata, ydata, p0=None):
    ((x0, y0, m1, m2), r) = scipy.optimize.curve_fit(piecewise_linear_f, xdata, ydata, p0=p0)
    return x0, y0, m1, m2, r


def shifted_rect_hyperbola_f(x, m, c):
    return m + c/x

def fit_shifted_rect_hyperbola(xdata, ydata, p0=None):
    ((m, c), r) = scipy.optimize.curve_fit(shifted_rect_hyperbola_f, xdata, ydata, p0=p0)
    return m, c, r
    
# def min_sep_energy_based_f(t, m, c,): # assumes no change in target final product purity
#                                       # and no significant resulting change in final water purity
#                                       # when titer changes
#     return m + c*-np.log(1-t)/t

#%% 4-coeff 2-param shifted rect hyperbola
# can use this when productivity=constant

def shifted_rect_hyperbola_two_param(x, y, a, b, c, d):
    return a + b/x + c/y + d/(x*y)

def shifted_rect_hyperbola_two_param_for_fit(xy, a, b, c, d):
    x, y = xy
    return a + b/x + c/y + d/(x*y)

def fit_shifted_rect_hyperbola_two_param(xydata, zdata, p0=None):
    # x, y = xydata
    ((a, b, c, d), r) = scipy.optimize.curve_fit(shifted_rect_hyperbola_two_param_for_fit, 
                                                  xydata, 
                                                  zdata, 
                                                  p0=p0,
                                                  xtol=1e-9,
                                                  ftol=1e-9,
                                                  maxfev=1000,)
    return a, b, c, d, r

# def shifted_rect_hyperbola_two_param(x, y, a, b, c, d):
#     return a + b/x + c/y 

# def shifted_rect_hyperbola_two_param_for_fit(xy, a, b, c, d):
#     x, y = xy
#     return a + b/x + c/y 

# def fit_shifted_rect_hyperbola_two_param(xydata, zdata, p0=None):
#     # x, y = xydata
#     ((a, b, c, d), r) = scipy.optimize.curve_fit(shifted_rect_hyperbola_two_param_for_fit, 
#                                                  xydata, 
#                                                  zdata, 
#                                                  p0=p0,
#                                                  xtol=1e-9,
#                                                  ftol=1e-9,
#                                                  maxfev=1000,)
#     return a, b, c, d, r


#%% 5-coeff 2-param shifted rect hyperbola
# !!! This is necessary if either productivity/time is involved
# If including productivity, change e*y to e*y/p
# If including time, change e*y to e*time

# new: e/productivity (that term is indep of y and t)
# when prod=constant, can use 4-coeff 2-param shifted rect hyperbola
def shifted_rect_hyperbola_5_coeff_two_param(x, y, z, a, b, c, d, e):
    return a + b/x + c/y + d/(x*y) + e/z

def shifted_rect_hyperbola_5_coeff_two_param_for_fit(xyz, a, b, c, d, e):
    x, y, z = xyz
    return a + b/x + c/y + d/(x*y) + e/z

def fit_shifted_rect_hyperbola_5_coeff_two_param(xydata, zdata, p0=None):
    # x, y = xydata
    ((a, b, c, d, e), r) = scipy.optimize.curve_fit(shifted_rect_hyperbola_5_coeff_two_param_for_fit, 
                                                 xydata, 
                                                 zdata, 
                                                 p0=p0,
                                                 xtol=1e-9,
                                                 ftol=1e-9,
                                                 maxfev=1000,)
    return a, b, c, d, e, r

# old: e*t
# def shifted_rect_hyperbola_5_coeff_two_param(x, y, a, b, c, d, e):
#     return a + b/x + c/y + d/(x*y) + e*y

# def shifted_rect_hyperbola_5_coeff_two_param_for_fit(xy, a, b, c, d, e):
#     x, y = xy
#     return a + b/x + c/y + d/(x*y) + e*y

# def fit_shifted_rect_hyperbola_5_coeff_two_param(xydata, zdata, p0=None):
#     # x, y = xydata
#     ((a, b, c, d, e), r) = scipy.optimize.curve_fit(shifted_rect_hyperbola_5_coeff_two_param_for_fit, 
#                                                  xydata, 
#                                                  zdata, 
#                                                  p0=p0,
#                                                  xtol=1e-9,
#                                                  ftol=1e-9,
#                                                  maxfev=1000,)
#     return a, b, c, d, e, r


# temp: e*y*t
def shifted_rect_hyperbola_5_coeff_two_param(x, y, a, b, c, d, e):
    return a + b/x + c/y + d/(x*y) + e*x*y

def shifted_rect_hyperbola_5_coeff_two_param_for_fit(xy, a, b, c, d, e):
    x, y = xy
    return a + b/x + c/y + d/(x*y) + e*x*y

def fit_shifted_rect_hyperbola_5_coeff_two_param(xydata, zdata, p0=None):
    # x, y = xydata
    ((a, b, c, d, e), r) = scipy.optimize.curve_fit(shifted_rect_hyperbola_5_coeff_two_param_for_fit, 
                                                 xydata, 
                                                 zdata, 
                                                 p0=p0,
                                                 xtol=1e-9,
                                                 ftol=1e-9,
                                                 maxfev=1000,)
    return a, b, c, d, e, r

#%% 6-coeff 2-param shifted rect hyperbola

def shifted_rect_hyperbola_two_param_6_coeff(x, y, a, b, c, d, e, f):
    return a + b/x + c/y + d/(x*y) + e*x + f*y

def shifted_rect_hyperbola_two_param_6_coeff_for_fit(xy, a, b, c, d, e, f):
    x, y = xy
    return a + b/x + c/y + d/(x*y) + e*x + f*y

def fit_shifted_rect_hyperbola_two_param_6_coeff(xydata, zdata, p0=None):
    # x, y = xydata
    ((a, b, c, d, e, f), r) = scipy.optimize.curve_fit(shifted_rect_hyperbola_two_param_6_coeff_for_fit, xydata, zdata, p0=p0)
    return a, b, c, d, e, f, r


#%%
def fit_min_sep_energy_based_f(xdata, ydata, p0=None):
    ((m,c), r) = scipy.optimize.curve_fit(min_sep_energy_based_f, xdata, ydata, p0=p0)
    return m, c, r

def min_sep_energy(mol_prod, mol_water,):# assumes no change in target final product purity
                                      # and no significant resulting change in final water purity
                                      # when titer changes
    return (  (mol_prod*0 + 0) # separated prod
            + (mol_water*0 + 0) # separated water
            -
              (mol_prod*np.log((mol_prod/(mol_prod+mol_water))) + mol_water*np.log((mol_water/(mol_prod+mol_water)))) # broth
            )

def min_sep_energy_based_f(t, m, c,): # assumes no change in target final product purity
                                      # and no significant resulting change in final water purity
                                      # when titer changes
    mol_prod = 1
    mol_water = mol_prod/t
    return m + c*min_sep_energy(mol_prod, mol_water,)

# def min_sep_energy_based_f(t, m, c, n, final_prod_purity=1.): # = product mol frac (mol_prod/mol_total)
#     water_mol_frac = 1. - t
#     final_water_mol_frac = 1. - final_prod_purity
#     if final_prod_purity < 1.:
#         return m + c*n*(
#                         (water_mol_frac*np.log(water_mol_frac) + 
#                         t*np.log(t))
#                         -
#                         (final_water_mol_frac*np.log(final_water_mol_frac) + 
#                         final_prod_purity*np.log(final_prod_purity))
#                         )
#     else:
#         return m + c*n*(
#                         (water_mol_frac*np.log(water_mol_frac) + 
#                         t*np.log(t))
#                         )
    
# def min_sep_energy_based_f_no_mol(t, m, c, final_prod_purity=1.): # = product mol frac (mol_prod/mol_total)
#     water_mol_frac = 1. - t
#     final_water_mol_frac = 1. - final_prod_purity
#     if final_prod_purity < 1.:
#         return m + c*(
#                         (final_water_mol_frac*np.log(final_water_mol_frac) + 
#                         final_prod_purity*np.log(final_prod_purity))
#                         -
#                         (water_mol_frac*np.log(water_mol_frac) + 
#                         t*np.log(t))
#                         )
#     else:
#         return m + c*-(
#                         (water_mol_frac*np.log(water_mol_frac) + 
#                         t*np.log(t))
#                         )
    
# def fit_min_sep_energy_based_f(tdata, indicatordata, ndata, y, final_prod_purity, p0=None):
#     def f(t, m, c, n, final_prod_purity):
#         return min_sep_energy_based_f(t, m, c, n, final_prod_purity)
#     ((m, c), r) = scipy.optimize.curve_fit(f, tdata, indicatordata, p0=p0)
#     return a, b, substrate_sugars_mol, r

# def min_sep_energy_based_f(t, a, b, substrate_sugars_mol, y): # t must be in mol/mol_total, not g/L
#     mol_prod = substrate_sugars_mol*y
#     mol_water = mol_prod/t - mol_prod
#     mol_total = mol_prod + mol_water
#     return (a - 
#             b*(
#                 mol_water*np.log(mol_water/mol_total) + 
#                 mol_prod*np.log(mol_prod/mol_total)
#                 )
#             )

# def fit_min_sep_energy_based_f(tdata, indicatordata, y, p0=None):
#     def f(t, a, b, substrate_sugars_mol):
#         return min_sep_energy_based_f(t=t, a=a, b=b, substrate_sugars_mol=substrate_sugars_mol, y=y)
#     ((a, b, substrate_sugars_mol), r) = scipy.optimize.curve_fit(f, tdata, indicatordata, p0=p0)
#     return a, b, substrate_sugars_mol, r


#%% Utils to get evenly distributed samples of feasible TY points

# yields given and returned in g/g but loaded to spec in %theoretical using theo_max_yield
# titers in g/L
def get_feasible_TY_samples(yields, titers, steps, MPSP_sim_f, theo_max_yield):
    titers_lb = titers[0]
    titers_reverse = np.flip(titers, 0)
    titers_ubs = []
    for y, i in zip(yields, range(len(yields))):
        titer_curr = None
        print(f'Looking at yield={y} ...')
        for t in titers_reverse:
            titer_curr = t
            try:
                MPSP = MPSP_sim_f(y/theo_max_yield, t)
                if not np.isnan(MPSP): 
                    print(f'yield={y}, titer={t} is feasible. MPSP = {MPSP}.')
                    titers_ubs.append(t)
                    break
            except:
                # print(f'{y}, {t} is infeasible.')
                pass
        if len(titers_ubs) < i+1:
            # print(titers_ub, i+1)
            breakpoint()
            raise RuntimeError(f'At yield {y}, no titer was identified with non-nan MPSP (lowest titer checked: {titer_curr}).')
            
    titer_samples = []
    yields_samples = np.linspace(yields[0], yields[-1], steps[0])
    yts = []
    for y, titers_ub in zip(yields_samples, titers_ubs):
        titer_samples.append(np.linspace(titers_lb, titers_ub, steps[1]))
        for t in titer_samples[-1]:
            yts.append((y, t))
    
    return yts