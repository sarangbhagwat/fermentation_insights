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
        shifted_rect_hyperbola_two_param,\
        plot_across_param
from biorefineries.TAL._general_utils import add_metrics_to_unit_groups,\
        TEA_breakdown as general_TEA_breakdown
from matplotlib import pyplot as plt
import contourplots
import os

from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap
import chaospy

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


#%% TAL TRY
# from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
#                                                             yields, titers, titers_mol_per_mol_total, productivities,\
#                                                             colors, CABBI_green_colormap, get_rounded_str,\
#                                                             R302, spec,\
#                                                             TAL_maximum_viable_market_range as market_range


# inflection_product_yields = [1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
#                             1-R302.regular_microbe_conversion]

#%% 3-HP TRY
from biorefineries.HP.analyses.fermentation.TRY_analysis_cornstover_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                            yields, titers, productivities,\
                                                            colors, CABBI_green_colormap, get_rounded_str,\
                                                            AA_market_range as market_range,\
                                                            R302, spec, get_AA_MPSP

inflection_product_yields = [1-R302.regular_biomass_conversion]

#%% Succinic TRY
# from biorefineries.succinic.analyses.TRY_analysis import results_metric_1, results_metric_2, results_metric_3,\
#                                                             yields, titers, productivities,\
#                                                             colors, CABBI_green_colormap,\
#                                                             R302, spec, get_product_MPSP
#                                                             # SA_price_range as market_range,\
# market_range = []
# inflection_product_yields = [1-0.186]

#%% 2,3-BDO TRY
# steps = 40
# yields = np.linspace(0.05, 0.95, steps)
# titers = np.linspace(5., 210., steps)
# productivities = np.array([1.00])

# market_range = []

# inflection_product_yields = [
#     # 1-0.055-0.02,
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
    
#%% Fit across titer at alternative yields
product_ID = 'HP'
plot_MPSPs = True
print_coeffs = True
plot_coefficients_vs_yield = True
# ms, cs = [], []
ms_fit_coeffs, cs_fit_coeffs = {}, {}
indicator_array = np.array(results_metric_1)
# indicator_array = np.array(results_metric_5/(yields*results_metric_6))
fit_or_solve = 'fit'
indices_for_solution = [3,8]

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
            m,c,r = fit_straight_line(titers_curr, titers_curr*MPSPs)
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
            ax.plot(titers_curr, [m + c/t for t in titers_curr], label='fit') # resulting MPSP vs titer curve
            plt.legend()
            
            ax = plt.subplot(413)
            ax.set_ylabel(f'Residuals [$/kg {product_ID}]')
            ax.set_xlabel(f'{product_ID} titer [g/L]')
            ax.scatter(titers_curr, MPSPs - np.array([m + c/t for t in titers_curr]))
            
            ax = plt.subplot(414)
            ax.set_ylabel('Standardized residuals [-]')
            ax.set_xlabel(f'{product_ID} titer [g/L]')
            ax.scatter(titers_curr, (MPSPs - np.array([m + c/t for t in titers_curr]))/MPSPs)
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
        ax.set_xlim(0)
        ax.set_ylim(0)
        ax = plt.subplot(212)
        ax.set_ylabel(f'c * yield[% theoretical]')
        ax.set_xlabel(f'{product_ID} yield[% theoretical]')
        ax.scatter(yields, cs*yields, label='orig') # from MPSP*titer vs titer straight line fit
        ax.plot(yields, piecewise_linear_f(yields, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2']),
                label='fit') # from (orig)*yield vs yield piecewise linear fit
        # ax.scatter(yields, cs, label='orig') # from MPSP*titer vs titer straight line fit
        # ax.plot(yields, piecewise_linear_f(yields, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2'])/yields,
        #         label='fit') # from (orig)*yield vs yield piecewise linear fit
        ax.set_xlim(0)
        ax.set_ylim(0)
        
        plt.legend()
        plt.show()
        
        # Fit shifted rect hyperbola to c values within a single yield regime 
        c_m, c_c, r = fit_shifted_rect_hyperbola(yields, cs, p0=[0.5,0.5])
        print(c_m, c_c)
        
def MPSP_y_t_f(y, t, ms_fit_coeffs, cs_fit_coeffs):
    y_inv = 1./y
    yield_index = np.where(yields==y)[0][0]
    titer_index = np.where(titers==t)[0][0]
    m = piecewise_linear_f(y, ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2']) * y_inv
    # m = ms[yield_index]
    c = piecewise_linear_f(y, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2']) * y_inv
    # c = cs[yield_index]
    if np.isnan(results_metric_1[0, titer_index, yield_index]): return np.nan
    return shifted_rect_hyperbola_f(t, m, c)

def mechanistic_MPSP_titer_f(t, a1, a2, t1, n, b):
    return a1 + a2*(t1/t)**n + b/t

#%% For solutions to system of eqns, sample separately for each yield regime

yield_upper_bound_index_for_eval = indices_inflection_product_yields[0]+1 # only evaluates for yield regimes that occur up to this yield index

def get_samples(yields, 
                titers, 
                indicator_array, 
                sampling_type='random',
                n_samples=4, 
                seed=None,
                ind_ind_for_MPSP_y=1, 
                ind_ind_for_MPSP_t=0, 
                low_bound_index_yield=0, 
                low_bound_index_titer=0,
                yield_upper_bound_index_for_eval=len(yields),
                repeating_indices_allowed=False):
    
    indices = []
    ys, ts, mpsps = [], [], []
    
    np.random.seed(seed)

    for i in [g for g in indices_inflection_product_yields if g<yield_upper_bound_index_for_eval]:
        # find random indices with non-nan results
        if sampling_type=='random':
            indices.append([])
            ys.append([])
            ts.append([])
            mpsps.append([])
            while len(indices[-1])<n_samples:
                ind = [np.random.randint(low_bound_index_yield, i), np.random.randint(low_bound_index_titer, len(titers)-1)]
                # ind = [np.random.randint(low_bound_index_yield, i), 15]
                indicatorvals = indicator_array[0][ind[ind_ind_for_MPSP_y], ind[ind_ind_for_MPSP_t]]
                if not np.any(np.isnan(indicatorvals)):
                    # print(indices)
                    if (len(indices[-1])==0) or\
                        repeating_indices_allowed or\
                        (not (np.any([indices[-1][k][0]==ind[0] or indices[-1][k][1]==ind[1] for k in range(len(indices[-1]))]))):
                        indices[-1].append(ind)
                        ys[-1].append(yields[ind[0]])
                        ts[-1].append(titers[ind[1]])
                        mpsps[-1].append(indicatorvals)
        
        elif sampling_type=='latin_hypercube':
            indic = np.array(list(zip(
                            np.int64(np.round(chaospy.Uniform(low_bound_index_yield, i).sample(n_samples, rule='latin_hypercube', seed=seed),0)),
                            np.int64(np.round(chaospy.Uniform(low_bound_index_titer, len(titers)-1).sample(n_samples, rule='random', seed=seed),0))
                          )))
            indices.append(indic)
            ys.append(yields[indic[:,0]])
            ts.append(titers[indic[:,1]])
            mpsps.append(indicator_array[0][indic[:,1], indic[:,0]])
        else:
            raise ValueError(f'{sampling_type}')
                    
        low_bound_index_yield = i+1

    return np.array(ys), np.array(ts), np.array(mpsps), np.array(indices)
    
ys, ts, mpsps, indices = get_samples(yields=yields[:yield_upper_bound_index_for_eval], titers=titers, indicator_array=indicator_array,
                                     yield_upper_bound_index_for_eval=yield_upper_bound_index_for_eval,
                                     repeating_indices_allowed=False)

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

#%% Solve using minimize (with constraints)

from scipy.optimize import minimize


def solve_for_coeffs(ys, ts, mpsps, yield_upper_bound_index_for_eval, print_results=True):
    avals, bvals, cvals, dvals = [], [], [], []
    
    for ydata, tdata, mpspdata in zip(ys, ts, mpsps):
        def minimize_obj_f(p):
            a, b, c, d = p
            # return np.abs(np.sum([
            #         a + b/ydata[0] + c/tdata[0] + d/(ydata[0]*tdata[0]) - mpspdata[0],
            #         a + b/ydata[1] + c/tdata[1] + d/(ydata[1]*tdata[1]) - mpspdata[1],
            #         a + b/ydata[2] + c/tdata[2] + d/(ydata[2]*tdata[2]) - mpspdata[2],
            #         a + b/ydata[3] + c/tdata[3] + d/(ydata[3]*tdata[3]) - mpspdata[3],
            #         ]))

            # return np.abs(a + b/ydata[0] + c/tdata[0] + d/(ydata[0]*tdata[0]) - mpspdata[0])\
            #      + np.abs(a + b/ydata[1] + c/tdata[1] + d/(ydata[1]*tdata[1]) - mpspdata[1])\
            #      + np.abs(a + b/ydata[2] + c/tdata[2] + d/(ydata[2]*tdata[2]) - mpspdata[2])\
            #      + np.abs(a + b/ydata[3] + c/tdata[3] + d/(ydata[3]*tdata[3]) - mpspdata[3])
        
            return np.abs(np.sum([
                    a + b/ydata[i] + c/tdata[i] + d/(ydata[i]*tdata[i]) - mpspdata[i]
                    for i in range(len(ydata))
                    ]))
            
            # return np.sum([
            #         np.abs(a + b/ydata[i] + c/tdata[i] + d/(ydata[i]*tdata[i]) - mpspdata[i])
            #         for i in range(len(ydata))
            #         ])
        cons = (
                # {'type': 'ineq', 'fun': lambda x:  x[0]},
    
                # {'type': 'ineq', 'fun': lambda x:  x[1]},
    
                # {'type': 'ineq', 'fun': lambda x:  x[2]},
                
                # {'type': 'ineq', 'fun': lambda x:  x[3]},
                
                {'type': 'eq', 'fun': lambda x: x[0] + x[1]/ydata[0] + x[2]/tdata[0] + x[3]/(ydata[0]*tdata[0]) - mpspdata[0]},
                {'type': 'eq', 'fun': lambda x: x[0] + x[1]/ydata[1] + x[2]/tdata[1] + x[3]/(ydata[1]*tdata[1]) - mpspdata[1]},
                # {'type': 'eq', 'fun': lambda x: x[0] + x[1]/ydata[2] + x[2]/tdata[2] + x[3]/(ydata[2]*tdata[2]) - mpspdata[2]},
                {'type': 'eq', 'fun': lambda x: x[0] + x[1]/ydata[3] + x[2]/tdata[3] + x[3]/(ydata[3]*tdata[3]) - mpspdata[3]},
                
                )
        
        res = minimize(minimize_obj_f, (1, 1, 1, 1), 
                       method='SLSQP', 
                       bounds=(
                           (0, None), # a
                           (0, None), # b
                           (0, None), # c
                           (0, None), # d
                           ),
                        constraints=cons,
                        tol=1e-8,
                        options={
                            'maxiter':1e3,
                            'disp':print_results,
                            'ftol':1e-8,
                            }
                       )
        a, b, c, d =  res.x
        if print_results:
            print(a, b, c, d)
            print(minimize_obj_f((a, b, c, d)))
            print(mpspdata)
            print(np.array([shifted_rect_hyperbola_two_param(yt[0], yt[1], a, b, c, d) for yt in zip(ydata, tdata)]))
            print('\n')
        try:
            rsq_check = get_Rsq(mpspdata, np.array([shifted_rect_hyperbola_two_param(yt[0], yt[1], a, b, c, d) for yt in zip(ydata, tdata)]))
            if print_results: print(rsq_check)
        except Exception as e:
            print('Could not compute R^2 value: {str(e)}.')
        avals.append(a)
        bvals.append(b)
        cvals.append(c)
        dvals.append(d)
    
    Rsq = (get_Rsq(indicator_array[:, :, :yield_upper_bound_index_for_eval], np.array([[[MPSP_f(y, t, avals, bvals, cvals, dvals) for y in yields[:yield_upper_bound_index_for_eval]] for t in titers]])))

    return Rsq, res.success, avals, bvals, cvals, dvals

Rsq, success, avals, bvals, cvals, dvals = solve_for_coeffs(ys=ys, ts=ts, mpsps=mpsps, yield_upper_bound_index_for_eval=yield_upper_bound_index_for_eval)

#%% Coefficient of determination

print(get_Rsq(indicator_array, np.array([[[MPSP_y_t_f(y,t,ms_fit_coeffs, cs_fit_coeffs) for y in yields] for t in titers]])))

#%%
print(get_Rsq(indicator_array[:, :, :yield_upper_bound_index_for_eval], np.array([[[MPSP_f(y, t, avals, bvals, cvals, dvals) for y in yields[:yield_upper_bound_index_for_eval]] for t in titers]])))

#%% Sampling-based uncertainty

def sample_and_solve(yields, 
                titers, 
                indicator_array, 
                sampling_type='random',
                n_samples=4, 
                seed=None,
                ind_ind_for_MPSP_y=1, 
                ind_ind_for_MPSP_t=0, 
                low_bound_index_yield=0, 
                low_bound_index_titer=0,
                yield_upper_bound_index_for_eval=indices_inflection_product_yields[0]+1,
                repeating_indices_allowed=False):
    ys, ts, mpsps, indices = get_samples(yields=yields[:yield_upper_bound_index_for_eval], titers=titers, indicator_array=indicator_array,
                                         yield_upper_bound_index_for_eval=yield_upper_bound_index_for_eval,
                                         repeating_indices_allowed=repeating_indices_allowed)
    # print(indices)
    # print(ys)
    Rsq, success, avals, bvals, cvals, dvals = solve_for_coeffs(ys=ys, ts=ts, mpsps=mpsps, yield_upper_bound_index_for_eval=yield_upper_bound_index_for_eval,
                                                       print_results=False)
    # print(Rsq, avals, bvals, cvals, dvals)
    return Rsq, success, avals, bvals, cvals, dvals
    
# For only 0th yield regime (no coproduct yield decrease)
unc_Rsqs = []
unc_as, unc_bs, unc_cs, unc_ds = [], [], [], []
n_solves = 10000
i=0
fail_ct = 0
while len(unc_Rsqs)<n_solves:
    Rsq, success, avals, bvals, cvals, dvals = sample_and_solve(yields=yields, 
                    titers=titers, 
                    indicator_array=indicator_array, 
                    sampling_type='random',
                    n_samples=4, 
                    seed=None,
                    ind_ind_for_MPSP_y=1, 
                    ind_ind_for_MPSP_t=0, 
                    low_bound_index_yield=0, 
                    low_bound_index_titer=0,
                    yield_upper_bound_index_for_eval=indices_inflection_product_yields[0]+1)
    unc_Rsqs.append(Rsq)
    unc_as.append(avals[0])
    unc_bs.append(bvals[0])
    unc_cs.append(cvals[0])
    unc_ds.append(dvals[0])
    if not success: fail_ct +=1
    i+=1
    if i%100==0:
        print(str(i)+f'/{n_solves}')
unc_Rsqs, unc_as, unc_bs, unc_cs, unc_ds = np.array(unc_Rsqs), np.array(unc_as), np.array(unc_bs), np.array(unc_cs), np.array(unc_ds)


plt.hist(unc_Rsqs[np.where(unc_Rsqs>0)], bins=100)
plt.show()
plt.hist(unc_as, bins=100)
plt.show()
plt.hist(unc_bs, bins=100)
plt.show()
plt.hist(unc_cs, bins=100)
plt.show()
plt.hist(unc_ds, bins=100)
plt.show()

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
