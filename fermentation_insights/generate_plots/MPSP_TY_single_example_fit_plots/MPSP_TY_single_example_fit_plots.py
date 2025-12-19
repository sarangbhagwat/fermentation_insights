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
import pandas as pd

from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
import chaospy
import itertools
import pickle
import dill
import imageio

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

product_IDs = [
                'TAL', 'TAL_SA', 
                'HP', 'HP_neutral', 
                'HP_hexanol', 'HP_neutral_hexanol', 
                'succinic', 'succinic_neutral',
                ]
feedstock_IDs = [
                 'glucose', 
                 'sugarcane', 
                 'corn', 
                 'cornstover',
                 ]

all_filenames = list(itertools.product(product_IDs, feedstock_IDs))

#%%
def format_ax(ax, x_ticks, y_ticks,):
    ax.tick_params( direction = 'inout' , which='both')
    ax.set_xlim((x_ticks[0], x_ticks[-1]))
    ax.set_ylim((y_ticks[0], y_ticks[-1]))
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    
    
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
    
#%% Load TRY

def get_MPSP_TY_single_example_fit_plots(
                    product, feedstock, additional_tag='', 
                    plot_MPSP_y_t=True, 
                    save_to='C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results'):
    os.chdir(save_to)
    product_ID = product
    # additional_tag = '0.5x_baselineprod'
    # additional_tag = 'neutral'
    # feedstock = 'cornstover'
    
    filename = None
    if additional_tag: 
        filename = f'{product}_{additional_tag}_{feedstock}'
    else:
        filename = f'{product}_{feedstock}'
    yields = np.load(f'{filename}_yields.npy')
    titers = np.load(f'{filename}_titers.npy')
    productivities = np.load(f'{filename}_productivities.npy')
    results_metric_1 = np.load(f'{filename}_MPSP.npy')
    results_metric_2 = np.load(f'{filename}_GWP.npy')
    results_metric_3 = np.load(f'{filename}_FEC.npy')
    inflection_product_yields = np.load(f'{filename}_inflection_product_yields.npy')
    recoveries = np.load(f'{filename}_recoveries.npy')
    
    
    # Recovery from after fermentation -> final product (incl. after catalytic upgrading)
    recovery = np.median(recoveries[np.where(~np.isnan(recoveries))].flatten())
    
    # Updated this to load recoveries over TRY and get median
    # recovery = None # recoveries in kg-final-product recovered / kg-fermentation-product-produced (where fermentation-product is the non-neutralized form)
    # if product=='HP': recovery = 0.6002
    # elif product=='HP_neutral': recovery = 0.6235
    # elif product=='HP_hexanol': recovery = 0.7070
    # elif product=='HP_neutral_hexanol': recovery = 0.7247
    # elif product=='TAL': recovery = 0.7191
    # elif product=='TAL_SA': recovery = 0.7078
    # elif product=='succinic': recovery = 0.9830
    
    # theoretical max fermentation yield, if loaded yields are in %theoretical max
    theo_max_yield = None # kg-fermentation-product-produced / kg-glucose-eq.
    if product in ('HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol'): theo_max_yield = 1.
    elif product in ('TAL', 'TAL_SA'): theo_max_yield = 0.4667
    elif product in ('succinic', 'succinic_neutral'): theo_max_yield = 1.311
    
    yields*=theo_max_yield
    inflection_product_yields*=theo_max_yield
    
    true_g = recovery

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
    
    #%%
    if not yields[-1] in inflection_product_yields: 
        inflection_product_yields.append(yields[-1])
        indices_inflection_product_yields.append(len(yields)-1)
          
    #%%
    indicator_array = np.array(results_metric_1)
    yield_upper_bound_index_for_eval = indices_inflection_product_yields[0]+1 # only evaluates for yield regimes that occur up to this yield index
    
    indicator_array_for_eval = indicator_array[:, :, :yield_upper_bound_index_for_eval]
    yields_for_eval = yields[:yield_upper_bound_index_for_eval]
    titers_for_eval = titers
    
    # np.save(f'{filename}_yields_for_eval.npy', yields_for_eval)
    # np.save(f'{filename}_titers_for_eval.npy', titers_for_eval)
    # np.save(f'{filename}_MPSP_for_eval.npy', indicator_array_for_eval)
    
    #%%
    def MPSP_f(y, t, avals, bvals, cvals, dvals): # use the coefficients solved for the corresponding yield regime
        for ii, i, iy in zip(range(len(indices_inflection_product_yields)),
                             indices_inflection_product_yields, 
                             inflection_product_yields):
            if y<=iy:
                return shifted_rect_hyperbola_two_param(y, t, avals[ii], bvals[ii], cvals[ii], dvals[ii],)
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
    
    
    #%% Across titers at alternative yields
    plot_MPSPs = True
    include_top_bar = True
    fontname={'fontname':'Arial Unicode'}
    z_ticks = yield_ticks = np.linspace(0., 0.8, 9)
    keep_frames = True
    print_coeffs = True
    plot_coefficients_vs_yield = True
    # ms, cs = [], []
    ms_fit_coeffs, cs_fit_coeffs = {}, {}
    # indicator_array = np.array(results_metric_5/(yields*results_metric_6))

    
    for i in range(len(indicator_array)):
        p = productivities[i]
        ms, cs = [], []
        t_M, t_ = [], []
        for j in range(len(indicator_array[0,0,:])):
            y = yields_for_eval[j]
            
            where_nan = np.where(np.isnan(indicator_array[0,:,j]))
            to_index = np.min(where_nan) if list(where_nan[0]) else -1 # lowest titer index at which MPSP=nan 
                                                                            # (i.e., titer is low enough at lower indices than this)
            MPSPs = indicator_array[0,:,j][:to_index]
            titers_curr = titers_for_eval[:to_index]
            
            m,c,r = None,None,None
            
            m,c,r = fit_straight_line(titers_curr, titers_curr*MPSPs)

                
            ms.append(m)
            cs.append(c)
            # if print_coeffs: print(y, ': ', m, c)
            t_M.append(titers_curr*MPSPs)
            t_.append(titers_curr)
            
        if plot_MPSPs:
            # fig, axs = plt.subplots(1, 1, constrained_layout=True,)
            
            for tmindex, T_, T_M in zip(range(len(t_)), t_, t_M):
                
                fig, axs = plt.subplots(2, 1, constrained_layout=True, 
                                        gridspec_kw={'height_ratios': [1, 20]})
                fig.set_figwidth(6)
                fig.set_figheight(6)
                if include_top_bar:
                    ax = axs[0]
                    a = [yields_for_eval[tmindex]]
                    ax.hlines(1,1,1)
                    ax.set_xlim(min(z_ticks), max(z_ticks))
                    ax.set_ylim(0.5,1.5)
                    ax.xaxis.tick_top()
                    ax.set_xlabel("Yield [g/g]",  
                                      # fontsize=12,
                                      **fontname)
                    ax.set_xticks(z_ticks,
                                  **fontname)
                
                    y = np.ones(np.shape(a))
                    ax.plot(a,y,
                            color='blue', 
                            marker='v',
                            ms = 7,)
                    ax.axes.get_yaxis().set_visible(False)
                    # ax.tick_params(labelsize=axis_tick_fontsize)
                    
                ax = axs[1]
                
                ax.set_xlim(0, 200)
                ax.set_ylim(0, 800)
                ax.set_ylabel(f'MPSP [$/kg {product_ID}] * titer[g/L]')
                ax.set_xlabel(f'{product_ID} titer [g/L]')
                ax.scatter(T_, T_M, label='simulated', color='blue', 
                           # alpha=tmindex/len(t_), 
                           marker='.',
                           )
                plt.savefig(f'MPSP_TY_single_example_frame_{tmindex}.png', 
                            transparent = False,  
                            facecolor = 'white',
                            bbox_inches='tight',
                            dpi=300,
                            )                                
                plt.close()
            frames = []
            for tmindex in range(len(t_)):
                image = imageio.v2.imread(f'MPSP_TY_single_example_frame_{tmindex}.png',)
                frames.append(image)
            
            
            # if n_loops==('inf' or 'infinite' or 'infinity' or np.inf):
            # breakpoint()
            imageio.mimsave('MPSP_TY_single_example_plot.gif',
                            frames,
                            fps=2,
                            loop=10,
                            ) 
                # frames.reverse()
                # imageio.mimsave('./' + 'reverse_'+animated_contourplot_filename + '.gif',
                #                 frames,
                #                 fps=fps,
                #                 ) 

            if not keep_frames:
                for tmindex in range(len(t_)):
                    os.remove(f'MPSP_TY_single_example_frame_{tmindex}.png',)
                 
    
            # ax.plot(titers_curr, [m*t + c for t in titers_curr], label='fit') # MPSP * titer vs titer straight line fit
                

                
                
        ms, cs = np.array(ms), np.array(cs)
        
        # Piecewise linear fit for m vs yield
        ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2'], ms_fit_coeffs['r'] = fit_piecewise_linear(yields, ms*yields, p0=[0.5,0.5,1.,1.])
        cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2'], cs_fit_coeffs['r'] = fit_piecewise_linear(yields, cs*yields, p0=[0.5,0.5,1.,1.])
        
        if plot_coefficients_vs_yield:
            fig = plt.figure()
            fig.set_figwidth(6)
            fig.set_figheight(6)
            # coeff*yield vs yield
            ax = plt.subplot(211)
            
            ax.set_xticks(z_ticks,
                          **fontname)
            ax.set_xlim(z_ticks[0], z_ticks[-1])
            # ax.set_ylim(0, 800)
            ax.set_ylabel(f'm * yield[% theoretical]')
            ax.set_xlabel(f'{product_ID} yield[% theoretical]')
            ax.scatter(yields[:-4], ms[:-4]*yields[:-4], label='orig',
            marker='.', color='blue') # from MPSP*titer vs titer straight line fit
            # ax.plot(yields, piecewise_linear_f(yields, ms_fit_coeffs['x0'], ms_fit_coeffs['y0'], ms_fit_coeffs['m1'], ms_fit_coeffs['m2']),
            #         label='fit') # from (orig)*yield vs yield piecewise linear fit
            
            ax = plt.subplot(212)
            ax.set_xticks(z_ticks,
                          **fontname)
            ax.set_xlim(z_ticks[0], z_ticks[-1])
            # ax.set_ylim(0, 800)
            ax.set_ylabel(f'c * yield[% theoretical]')
            ax.set_xlabel(f'{product_ID} yield[% theoretical]')
            ax.scatter(yields[:-4], cs[:-4]*yields[:-4], label='orig',
            marker='.', color='blue') # from MPSP*titer vs titer straight line fit
            # ax.plot(yields, piecewise_linear_f(yields, cs_fit_coeffs['x0'], cs_fit_coeffs['y0'], cs_fit_coeffs['m1'], cs_fit_coeffs['m2']),
            #         label='fit') # from (orig)*yield vs yield piecewise linear fit
            
            # plt.legend()
            # plt.show()
            plt.savefig(f'slope_intercept_MPSP_TY_single_example.png', 
                        transparent = False,  
                        facecolor = 'white',
                        bbox_inches='tight',
                        dpi=300,
                        )  
            plt.close()
