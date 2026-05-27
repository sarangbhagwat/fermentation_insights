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
import pandas as pd

from openpyxl import load_workbook
from pathlib import Path

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

#%% Multivariate linear regression

def multivariate_linear_eq(x, y, a, b, c):
    return a + b*x + c*y

def multivariate_linear_eq_for_fit(xy, a, b, c):
    x, y = xy
    return a + b*x + c*y
    
def fit_multivariate_linear_eq(xydata, zdata, p0=None):
    # x, y = xydata
    ((a, b, c), r) = scipy.optimize.curve_fit(multivariate_linear_eq_for_fit, xydata, zdata, p0=p0)
    return a, b, c, r

#%% Alternative regression suggested by a reviewer

def alternative_eq(x, y, a, b, c, d):
    return a + b*x + c*y + d*x*y

def alternative_eq_for_fit(xy, a, b, c, d):
    x, y = xy
    return a + b*x + c*y + d*x*y
    
def fit_alternative_eq(xydata, zdata, p0=None):
    # x, y = xydata
    ((a, b, c, d), r) = scipy.optimize.curve_fit(alternative_eq_for_fit, xydata, zdata, p0=p0)
    return a, b, c, d, r

#%% Utils to get evenly distributed samples of feasible TY points

# yields given and returned in g/g but loaded to spec in %theoretical using theo_max_yield
# titers in g/L
def get_feasible_TY_samples(yields, titers, steps, MPSP_sim_f, theo_max_yield):
    titers_lb = titers[0]
    titers_reverse = np.flip(titers, 0)
    titers_ubs = []
    
    yields_samples = np.linspace(yields[0], yields[-1], steps[0])
    
    for y, i in zip(yields_samples, range(len(yields_samples))):
        # titers_ubs.append([])
        titer_curr = None
        print(f'Looking at yield={y} ...')
        for t in titers_reverse:
            titer_curr = t
            try:
                MPSP = MPSP_sim_f(y/theo_max_yield, t)
                if not np.isnan(MPSP): 
                    print(f'For yield={y}, titer={t} is the highest titer feasible within the given range. MPSP = {MPSP}.')
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
    yts = []
    for y, titers_ub in zip(yields_samples, titers_ubs):
        titer_samples.append(np.linspace(titers_lb, titers_ub, steps[1]))
        for t in titer_samples[-1]:
            yts.append((y, t))
    
    return yts

#%% Temporary functions to get all pairings from HXN

def get_pair_indices_from_unit_ID(ID):
    k = ID.replace('HX_', '').replace('_hs', '').replace('_cs', '')
    i1, i2 = int(k[:k.index('_')]), int(k[k.index('_')+1:])
    return i1, i2

def get_pairings(HXN, owner=False):
    hxs = HXN.new_HXs
    orig_hus = HXN.original_heat_utils
    pairings = []
    for hx in hxs:
        i1, i2 = get_pair_indices_from_unit_ID(hx.ID)
        if not owner:
            pairings.append((orig_hus[i1].unit, orig_hus[i2].unit))
        else:
            pairings.append((orig_hus[i1].unit.owner, orig_hus[i2].unit.owner))
    return pairings

#%%

def save_xyz_grid_to_excel(x_values, y_values, z_values, output_path, sheet_name, x_title, y_title, z_title):
    """
    Save x, y, z grid data to an Excel sheet.

    Layout:
        - First row contains x-axis values
        - First column contains y-axis values
        - Middle cells contain z-values

    Parameters
    ----------
    x_values : 1D array-like
        X-axis values, length = number of columns in z_values

    y_values : 1D array-like
        Y-axis values, length = number of rows in z_values

    z_values : 2D array-like
        Z-values with shape (len(y_values), len(x_values))

    output_path : str
        Path to output Excel file, e.g. "output.xlsx"
    
    sheet_name : str
        Name of the Excel sheet to which to write.
    
    x_titles : str
        X-axis title.

    y_titles : str
        Y-axis title.

    z_titles : str
        Z-axis title.
    """

    x_values = np.asarray(x_values)
    y_values = np.asarray(y_values)
    z_values = np.asarray(z_values)

    if z_values.shape != (len(y_values), len(x_values)):
        raise ValueError(
            f"z_values must have shape ({len(y_values)}, {len(x_values)}), "
            f"but got {z_values.shape}"
        )

    # Create DataFrame with y-values as row index and x-values as column headers
    df = pd.DataFrame(
        data=z_values,
        index=y_values,
        columns=x_values
    )

    # Optional: name the index/column corner cell
    df.index.name = f"{y_title} \\ {x_title}"

    # Save to Excel
    
    with pd.ExcelWriter(output_path,
        engine="openpyxl",
        mode="a",
        if_sheet_exists="overlay", # or "new", "replace"
        ) as writer:
        
        df.to_excel(
            writer,
            sheet_name=sheet_name,
            startrow=1
        )

        # Write z axis title to cell A1
        worksheet = writer.sheets[sheet_name]
        worksheet["A1"] = z_title


#%%
def sort_excel_sheets_alphabetically(input_path, output_path=None):
    input_path = Path(input_path)

    if output_path is None:
        output_path = input_path.with_stem(input_path.stem + "_sorted")

    wb = load_workbook(input_path)

    wb._sheets.sort(key=lambda ws: ws.title.casefold())

    wb.save(output_path)
    return output_path
