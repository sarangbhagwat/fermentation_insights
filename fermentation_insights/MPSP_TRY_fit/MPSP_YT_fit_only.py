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

from labellines import labelLine, labelLines

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

def get_all_MPSP_yt_fit():
    not_found = []
    for i in all_filenames:
        print(i[0], i[1])
        for additional_tag in ('', '0.2bp', '5.0bp'):
            try:
                get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag=additional_tag, plot_MPSP_y_t=False)
            except:
                not_found.append((i[1], i[0], additional_tag))
        

        if i[0]+'_'+i[1] in ['TAL_SA_sugarcane']:
            get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag='0.2bp', plot_MPSP_y_t=False)
            # get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag='1.0bp', plot_MPSP_y_t=False)
            get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag='1.8bp', plot_MPSP_y_t=False)
            get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag='2.6bp', plot_MPSP_y_t=False)
            get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag='3.4bp', plot_MPSP_y_t=False)
            get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag='4.2bp', plot_MPSP_y_t=False)
            get_MPSP_yt_fit(product=i[0], feedstock=i[1], additional_tag='5.0bp', plot_MPSP_y_t=False)
    
    return not_found

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

def get_MPSP_yt_fit(product, feedstock, additional_tag='', 
                    plot_MPSP_y_t=False, 
                    external_data_fit_and_plot=False,
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
    
    #%% TRY simulations
    # #%% TAL TRY
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'TAL'
    # if run_simulations:
    #     from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_cornstover import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
    #                                                                 yields, titers, titers_mol_per_mol_total, productivities,\
    #                                                                 colors, CABBI_green_colormap, get_rounded_str,\
    #                                                                 R302, spec,\
    #                                                                 TAL_maximum_viable_market_range as market_range
        
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
    #                                           1-R302.regular_microbe_conversion])
    #     np.save('TAL_cornstover_MPSP', results_metric_1)
    #     np.save('TAL_cornstover_GWP', results_metric_2)
    #     np.save('TAL_cornstover_FEC', results_metric_3)
    #     np.save('TAL_cornstover_AOC', results_metric_5)
    #     np.save('TAL_cornstover_TCI', results_metric_6)
    #     np.save('TAL_cornstover_recoveries', results_metric_4)
    #     np.save('TAL_cornstover_inflection_product_yields', inflection_product_yields)
    #     np.save('TAL_cornstover_yields', yields)
    #     np.save('TAL_cornstover_titers', titers)
    #     np.save('TAL_cornstover_productivities', productivities)
    # else:
    #     results_metric_1 = np.load('TAL_cornstover_MPSP')
    #     results_metric_2 = np.load('TAL_cornstover_GWP')
    #     results_metric_3 = np.load('TAL_cornstover_FEC')
    #     results_metric_5 = np.load('TAL_cornstover_AOC')
    #     results_metric_6 = np.load('TAL_cornstover_TCI')
    #     inflection_product_yields = np.load('TAL_cornstover_inflection_product_yields')
                
    # #%% TAL TRY - SA
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'TAL-SA'
    # if run_simulations:
    #     from biorefineries.TAL.analyses.fermentation.TRY_analysis_SA_sugarcane import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
    #                                                                 yields, titers, titers_mol_per_mol_total, productivities,\
    #                                                                 colors, CABBI_green_colormap, get_rounded_str,\
    #                                                                 R302, spec
                                                                    
        
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
    #                                           1-R302.regular_microbe_conversion])
    #     np.save('TAL_SA_sugarcane_MPSP', results_metric_1)
    #     np.save('TAL_SA_sugarcane_GWP', results_metric_2)
    #     np.save('TAL_SA_sugarcane_FEC', results_metric_3)
    #     np.save('TAL_SA_sugarcane_AOC', results_metric_5)
    #     np.save('TAL_SA_sugarcane_TCI', results_metric_6)
    #     np.save('TAL_SA_sugarcane_recoveries', results_metric_4)
    #     np.save('TAL_SA_sugarcane_inflection_product_yields', inflection_product_yields)
    #     np.save('TAL_SA_sugarcane_yields', yields)
    #     np.save('TAL_SA_sugarcane_titers', titers)
    #     np.save('TAL_SA_sugarcane_productivities', productivities)
        
    # else:
    #     results_metric_1 = np.load('TAL_SA_sugarcane_MPSP')
    #     results_metric_2 = np.load('TAL_SA_sugarcane_GWP')
    #     results_metric_3 = np.load('TAL_SA_sugarcane_FEC')
    #     results_metric_5 = np.load('TAL_SA_sugarcane_AOC')
    #     results_metric_6 = np.load('TAL_SA_sugarcane_TCI')
    #     inflection_product_yields = np.load('TAL_SA_sugarcane_inflection_product_yields')
                
    
    # #%% 3-HP TRY
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'HP'
    # if run_simulations:
    #     from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
    #                                                                 yields, titers, productivities,\
    #                                                                 colors, CABBI_green_colormap, get_rounded_str,\
    #                                                                 AA_market_range as market_range,\
    #                                                                 R302, spec, get_AA_MPSP
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
        
    #     np.save('HP_glucose_MPSP', results_metric_1)
    #     np.save('HP_glucose_GWP', results_metric_2)
    #     np.save('HP_glucose_FEC', results_metric_3)
    #     np.save('HP_glucose_AOC', results_metric_4)
    #     np.save('HP_glucose_TCI', results_metric_5)
    #     np.save('HP_glucose_recoveries', results_metric_6)
    #     np.save('HP_glucose_inflection_product_yields', inflection_product_yields)
    #     np.save('HP_glucose_yields', yields)
    #     np.save('HP_glucose_titers', titers)
    #     np.save('HP_glucose_productivities', productivities)
    # else:
    #     results_metric_1 = np.load('HP_glucose_MPSP')
    #     results_metric_2 = np.load('HP_glucose_GWP')
    #     results_metric_3 = np.load('HP_glucose_FEC')
    #     results_metric_4 = np.load('HP_glucose_AOC')
    #     results_metric_5 = np.load('HP_glucose_TCI')
    #     inflection_product_yields = np.load('HP_glucose_inflection_product_yields')
    
    # #%% 3-HP Neutral TRY
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'HP'
    # if run_simulations:
    #     from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
    #                                                                 yields, titers, productivities,\
    #                                                                 colors, CABBI_green_colormap, get_rounded_str,\
    #                                                                 AA_market_range as market_range,\
    #                                                                 R302, spec, get_AA_MPSP
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
        
    #     np.save('HP_neutral_glucose_MPSP', results_metric_1)
    #     np.save('HP_neutral_glucose_GWP', results_metric_2)
    #     np.save('HP_neutral_glucose_FEC', results_metric_3)
    #     np.save('HP_neutral_glucose_AOC', results_metric_4)
    #     np.save('HP_neutral_glucose_TCI', results_metric_5)
    #     np.save('HP_neutral_glucose_recoveries', results_metric_6)
    #     np.save('HP_neutral_glucose_inflection_product_yields', inflection_product_yields)
    #     np.save('HP_neutral_glucose_yields', yields)
    #     np.save('HP_neutral_glucose_titers', titers)
    #     np.save('HP_neutral_glucose_productivities', productivities)
    # else:
    #     results_metric_1 = np.load('HP_neutral_glucose_MPSP')
    #     results_metric_2 = np.load('HP_neutral_glucose_GWP')
    #     results_metric_3 = np.load('HP_neutral_glucose_FEC')
    #     results_metric_4 = np.load('HP_neutral_glucose_AOC')
    #     results_metric_5 = np.load('HP_neutral_glucose_TCI')
    #     inflection_product_yields = np.load('HP_neutral_glucose_inflection_product_yields')
    
    # #%% 3-HP Hexanol TRY
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'HP'
    # if run_simulations:
    #     from biorefineries.HP.analyses.fermentation.TRY_analysis_corn_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
    #                                                                 yields, titers, productivities,\
    #                                                                 colors, CABBI_green_colormap, get_rounded_str,\
    #                                                                 AA_market_range as market_range,\
    #                                                                 R302, spec, get_AA_MPSP
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
        
    #     np.save('HP_hexanol_corn_MPSP', results_metric_1)
    #     np.save('HP_hexanol_corn_GWP', results_metric_2)
    #     np.save('HP_hexanol_corn_FEC', results_metric_3)
    #     np.save('HP_hexanol_corn_AOC', results_metric_4)
    #     np.save('HP_hexanol_corn_TCI', results_metric_5)
    #     np.save('HP_hexanol_corn_recoveries', results_metric_6)
    #     np.save('HP_hexanol_corn_inflection_product_yields', inflection_product_yields)
    #     np.save('HP_hexanol_corn_yields', yields)
    #     np.save('HP_hexanol_corn_titers', titers)
    #     np.save('HP_hexanol_corn_productivities', productivities)
    # else:
    #     results_metric_1 = np.load('HP_hexanol_corn_MPSP')
    #     results_metric_2 = np.load('HP_hexanol_corn_GWP')
    #     results_metric_3 = np.load('HP_hexanol_corn_FEC')
    #     results_metric_4 = np.load('HP_hexanol_corn_AOC')
    #     results_metric_5 = np.load('HP_hexanol_corn_TCI')
    #     inflection_product_yields = np.load('HP_hexanol_corn_inflection_product_yields')
    
    # #%% 3-HP Hexanol Neutral TRY
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'HP'
    # if run_simulations:
    #     from biorefineries.HP.analyses.fermentation.TRY_analysis_cornstover_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
    #                                                                 yields, titers, productivities,\
    #                                                                 colors, CABBI_green_colormap, get_rounded_str,\
    #                                                                 AA_market_range as market_range,\
    #                                                                 R302, spec, get_AA_MPSP
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
        
    #     np.save('HP_neutral_hexanol_cornstover_MPSP', results_metric_1)
    #     np.save('HP_neutral_hexanol_cornstover_GWP', results_metric_2)
    #     np.save('HP_neutral_hexanol_cornstover_FEC', results_metric_3)
    #     np.save('HP_neutral_hexanol_cornstover_AOC', results_metric_4)
    #     np.save('HP_neutral_hexanol_cornstover_TCI', results_metric_5)
    #     np.save('HP_neutral_hexanol_cornstover_recoveries', results_metric_6)
    #     np.save('HP_neutral_hexanol_cornstover_inflection_product_yields', inflection_product_yields)
    #     np.save('HP_neutral_hexanol_cornstover_yields', yields)
    #     np.save('HP_neutral_hexanol_cornstover_titers', titers)
    #     np.save('HP_neutral_hexanol_cornstover_productivities', productivities)
    # else:
    #     results_metric_1 = np.load('HP_neutral_hexanol_cornstover_MPSP')
    #     results_metric_2 = np.load('HP_neutral_hexanol_cornstover_GWP')
    #     results_metric_3 = np.load('HP_neutral_hexanol_cornstover_FEC')
    #     results_metric_4 = np.load('HP_neutral_hexanol_cornstover_AOC')
    #     results_metric_5 = np.load('HP_neutral_hexanol_cornstover_TCI')
    #     inflection_product_yields = np.load('HP_neutral_hexanol_cornstover_inflection_product_yields')
    
    # #%% Succinic TRY
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'succinic'
    # if run_simulations:
    #     from biorefineries.succinic.analyses.TRY_analysis_cornstover import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
    #                                                                 yields, titers, productivities,\
    #                                                                 colors, CABBI_green_colormap,\
    #                                                                 R302, spec, get_product_MPSP
    #                                                                 # SA_price_range as market_range,\
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     market_range = []
    #     inflection_product_yields = np.array([1-0.186])
    #     np.save('succinic_cornstover_MPSP', results_metric_1)
    #     np.save('succinic_cornstover_GWP', results_metric_2)
    #     np.save('succinic_cornstover_FEC', results_metric_3)
    #     # np.save('succinic_cornstover_AOC', results_metric_4)
    #     # np.save('succinic_cornstover_TCI', results_metric_5)
    #     np.save('succinic_cornstover_recoveries', results_metric_4)
    #     np.save('succinic_cornstover_inflection_product_yields', inflection_product_yields)
    #     np.save('succinic_cornstover_yields', yields)
    #     np.save('succinic_cornstover_titers', titers)
    #     np.save('succinic_cornstover_productivities', productivities)
    # else:
    #     results_metric_1 = np.load('succinic_cornstover_MPSP')
    #     results_metric_2 = np.load('succinic_cornstover_GWP')
    #     results_metric_3 = np.load('succinic_cornstover_FEC')
    #     recoveries = np.load('succinic_cornstover_recoveries')
    #     # results_metric_4 = np.load('succinic_cornstover_AOC')
    #     # results_metric_5 = np.load('succinic_cornstover_TCI')
    #     inflection_product_yields = np.load('succinic_cornstover_inflection_product_yields')
        
    # #%% Succinic Neutral TRY
    # run_simulations = True
    # import numpy as np
    # import os
    # product_ID = 'succinic'
    # if run_simulations:
    #     from biorefineries.succinic.analyses.TRY_analysis_glucose import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
    #                                                                 yields, titers, productivities,\
    #                                                                 colors, CABBI_green_colormap,\
    #                                                                 R302, spec, get_product_MPSP
    #                                                                 # SA_price_range as market_range,\
    #     os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    #     market_range = []
    #     inflection_product_yields = np.array([1-0.186])
    #     np.save('succinic_neutral_glucose_MPSP', results_metric_1)
    #     np.save('succinic_neutral_glucose_GWP', results_metric_2)
    #     np.save('succinic_neutral_glucose_FEC', results_metric_3)
    #     # np.save('succinic_neutral_glucose_AOC', results_metric_4)
    #     # np.save('succinic_neutral_glucose_TCI', results_metric_5)
    #     np.save('succinic_neutral_glucose_recoveries', results_metric_4)
    #     np.save('succinic_neutral_glucose_inflection_product_yields', inflection_product_yields)
    #     np.save('succinic_neutral_glucose_yields', yields)
    #     np.save('succinic_neutral_glucose_titers', titers)
    #     np.save('succinic_neutral_glucose_productivities', productivities)
    # else:
    #     results_metric_1 = np.load('succinic_neutral_glucose_MPSP')
    #     results_metric_2 = np.load('succinic_neutral_glucose_GWP')
    #     results_metric_3 = np.load('succinic_neutral_glucose_FEC')
    #     recoveries = np.load('succinic_neutral_glucose_recoveries')
    #     # results_metric_4 = np.load('succinic_neutral_glucose_AOC')
    #     # results_metric_5 = np.load('succinic_neutral_glucose_TCI')
    #     inflection_product_yields = np.load('succinic_neutral_glucose_inflection_product_yields')
        
        
    # #%% 2,3-BDO TRY
    # # product_ID = 'BDO'
    # # ## yield-titer indices should be flipped for BDO, as below
    # # product_ID = 'BDO'
    # # ind_ind_for_MPSP_y = 0
    # # ind_ind_for_MPSP_t = 1
    # # ##
    
    # # steps = 40
    # # yields = np.linspace(0.05, 0.95, steps)
    # # titers = np.linspace(5., 210., steps)
    # # productivities = np.array([1.00])
    
    # # market_range = []
    
    # # inflection_product_yields = [
    # #     1-0.055-0.02,
    # #     # 1-0.055,
    # #     ]
    
    # # os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results//BDO')
    
    # # arr = np.load('BDO_TRY_2021.11.9.18.57'+'.npy')
    
    # # # results_metric_1 = arr[:, :, 0, :].transpose()/907.185
    # # # results_metric_2 = arr[:, :, 1, :].transpose()
    # # # results_metric_3 = arr[:, :, 2, :].transpose()
    
    # # results_metric_1 = arr[:, :, 0, :]/907.185
    # # results_metric_2 = arr[:, :, 1, :]
    # # results_metric_3 = arr[:, :, 2, :]
    
    # # results_metric_1_new, results_metric_2_new, results_metric_3_new = [], [], []
    # # for i, j, k in zip(results_metric_1, results_metric_2, results_metric_3):
    # #     results_metric_1_new.append(i.flatten())
    # #     results_metric_2_new.append(j.flatten())
    # #     results_metric_3_new.append(k.flatten())
        
    # # results_metric_1 = np.array([results_metric_1_new])
    # # results_metric_2 = np.array([results_metric_2_new])
    # # results_metric_3 = np.array([results_metric_3_new])
    
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
    titers_mol_mol = False
    titers_gpL = None
    if titers_mol_mol:
        titers_gpL = titers.copy()
        titers = np.array(titers_mol_per_mol_total[0])
        
    #%%
    indicator_array = np.array(results_metric_1)
    yield_upper_bound_index_for_eval = indices_inflection_product_yields[0]+1 # only evaluates for yield regimes that occur up to this yield index
    
    indicator_array_for_eval = indicator_array[:, :, :yield_upper_bound_index_for_eval]
    yields_for_eval = yields[:yield_upper_bound_index_for_eval]
    titers_for_eval = titers
    
    np.save(f'{filename}_yields_for_eval.npy', yields_for_eval)
    np.save(f'{filename}_titers_for_eval.npy', titers_for_eval)
    np.save(f'{filename}_MPSP_for_eval.npy', indicator_array_for_eval)
    
    
    np.save(f'{filename}_recoveries_for_eval.npy', recoveries[:, :, :yield_upper_bound_index_for_eval])
    
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
                                                                               [1,1,1,1],
                                                                               # [0,0,0,0,0],
                                                                               )
    
    # print('\n')
    # print(filename)
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
            # g, aprime, bprime, cprime, dprime, eprime = x
            return sum([(x[1] +
                          x[2]/ys[i] + 
                          x[3]/ts[i] + 
                          x[4]/(ys[i]*ts[i]) - x[0]*mpsps[i])**2
                          for i in range(len(ys))])
            # return sum([(x[1]*y +
            #              x[2] + 
            #              x[3]*y/t + 
            #              x[4]*t - x[0]*mpsp*y)**2
            #              for y,t,mpsp in zip(ys,ts,mpsps)])
        
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
                        bounds=(
                            (100*216/227, 100*216/227), # g
                            (None, None), # aprime
                            (None, None), # bprime
                            (None, None), # cprime
                            (None, None), # dprime
                            ),
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
            print(aprime, bprime, cprime, dprime, g)
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
                                                                          [shifted_rect_hyperbola_two_param(yt[0], yt[1], a_fit, b_fit, c_fit, d_fit,) 
                                                                           for yt in zip(ys_fit[:2], ts_fit[:2])],
                                                                          a_fit, b_fit, c_fit, d_fit)
    
    
    #%% Relative impact
    
    import sympy as sp
    print('\n\n\n')
    y, t = sp.symbols('y t')
    a,b,c,d = a_fit, b_fit, c_fit, d_fit
    
    f = a + b*100/y +c/t + d*100/(y*t) # if yield increment to test is 0.01 g/g rather than 1 g/g, multiple b and d by 100
    
    # Partial derivative with respect to x
    df_dy = sp.diff(f, y) 
    print(df_dy)
    
    # Partial derivative with respect to y
    df_dt = sp.diff(f, t) 
    print(df_dt)
    
    rel_impact = df_dy/df_dt
    print(rel_impact)
    
    # rel_impact_f = sp.utilities.lambdify((y,t), rel_impact)
    
    # lambda function without using sympy
    def rel_impact_f(y, t, y_increment=0.01, t_increment=1):
        # if yield increment to test is 0.01 g/g rather than 1 g/g, multiple b and d by 100
        # if titer increment to test is 5 g/L rather than 1 g/L, multiple c and d by 0.2
        # yes, if changing both increments, multiply d by both
        # expects y in units of (g/g)/(y_increment), and t in units of (g/L)/(t_increment)
        b_ri = b/y_increment
        c_ri = c/t_increment
        d_ri = d/(y_increment*t_increment)
        return ((-b_ri/y**2 - d_ri/(t*(y**2)))/
                (-c_ri/t**2 - d_ri/((t**2)*y)))
    
    with open(filename+'_RI_f_MPSP.pkl', 'wb') as f:
        dill.dump(rel_impact_f, f)
        # dill.dumps(rel_impact_f, recurse=True)
    
    with open(filename+'_RI_f_MPSP.pkl', 'rb') as f:
        rel_impact_f = dill.load(f)
    
    y_increment=0.01 
    t_increment=1
    
    ys_ri = yields_for_eval/y_increment
    ts_ri = titers_for_eval/t_increment

    rel_impact_arr = np.array([[[rel_impact_f(yi, ti) for yi in ys_ri] for ti in ts_ri]])
    
    rel_impact_arr[np.where(np.isnan(indicator_array_for_eval))] = np.nan
    
    np.save(f'{filename}_RI_MPSP.npy', rel_impact_arr)
    np.save(f'{filename}_yields_RI.npy', ys_ri)
    np.save(f'{filename}_titers_RI.npy', ts_ri)
    
    #%% Save coefficients
    
    np.save(
            filename+'_coefficients',
            np.array(
                    [a_fit*true_g, 
                     b_fit*true_g, 
                     c_fit*true_g, 
                     d_fit*true_g,
                     true_g,
                     get_Rsq(indicator_array_for_eval, 
                                  np.array([[[MPSP_f(y, t, [a_fit], [b_fit], [c_fit], [d_fit]) 
                                              for y in yields_for_eval] 
                                                 for t in titers_for_eval]]))
                    ]
                    ),
                    )
    
    #%% Plot stuff
    yields_for_plot = yields_for_eval
    titers_for_plot = titers_for_eval
    indicator_array_for_plot = indicator_array_for_eval.tolist()
    
    #%% MPSP_fit
    fit_indicator_array_for_plot = [[[MPSP_f(y, t, [a_fit], [b_fit], [c_fit], [d_fit]) 
                                   for y in yields_for_plot] 
                                   for t in titers_for_plot]]
    
    # for i in range(len(indicator_array)):
    #     for j in range(len(indicator_array[i])):
    #         for k in range(len(indicator_array[i][j])):
    #             # print(k)
    #             if k>=(len(indicator_array_for_plot[i][j])):
    #                 indicator_array_for_plot[i][j].append(np.nan)
    #             if k>=(len(fit_indicator_array_for_plot[i][j])):
    #                 fit_indicator_array_for_plot[i][j].append(np.nan)
    
    
    indicator_array_for_plot = np.array(indicator_array_for_plot)
    fit_indicator_array_for_plot = np.array(fit_indicator_array_for_plot)
    fit_indicator_array_for_plot[np.where(np.isnan(indicator_array_for_plot))] = np.nan
    
    np.save(f'{filename}_MPSP_fit.npy', fit_indicator_array_for_plot)
    
    
    #%% Contribution % of d term
    def dcon_magnitude_f(y, t, d):
        return d/(y*t)
    dcon = [[[ (        dcon_magnitude_f(y, t, d_fit)/
                        MPSP_f(y, t, [a_fit], [b_fit], [c_fit], [d_fit])
                )
                                   for y in yields_for_plot] 
                                   for t in titers_for_plot]]
    dcon = np.array(dcon)
    
    dcon[np.where(np.isnan(indicator_array_for_plot))] = np.nan
    
    dcon_non_nan = dcon[np.where(~np.isnan(dcon))]
    
    print('D and DCON min and max:', d_fit, np.min(dcon_non_nan), np.max(dcon_non_nan))
    
    np.save(f'{filename}_dcon.npy', dcon)
    
    #%% Plots
    if plot_MPSP_y_t:
        # Parameters analyzed across
    
        x_label = r"$\bfYield$" # title of the x axis
        x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
        # x_ticks = [20, 30, 40, 50, 60, 70, 80, 90, 100]
        
        round_up_yield_ubound = True if np.round(yields_for_plot[-1],1)<yields_for_plot[-1] else False
        round_down_yield_lbound = True if np.round(yields_for_plot[0],1)>yields_for_plot[0] else False
        
        x_ticks = np.arange(np.round(yields_for_plot[0],1)-0.1*round_down_yield_lbound, 
                            np.round(yields_for_plot[-1],1)+0.1*round_up_yield_ubound+1e-5, 0.1)[::1]
        
        y_label = r"$\bfTiter$" # title of the y axis
        y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
        
        round_up_titer_ubound = True if np.round(titers_for_plot[-1],-1)<titers_for_plot[-1] else False
        round_down_titer_lbound = True if np.round(titers_for_plot[0],-1)>titers_for_plot[0] else False
        
        y_ticks = np.arange(np.round(titers_for_plot[0],-1)-10*round_down_titer_lbound, 
                            np.round(titers_for_plot[-1],-1)+10*round_up_titer_ubound+1e-5, 10)[::2]
    
    
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
        
        rel_impact_w_label = r"$\bfRelative$" + " " + r"$\bfImpact$"
        rel_impact_units = r"$\mathrm{kg-sugars}\cdot\mathrm{L-broth}^{-1}$"
        
        fps = 3
        axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
        default_fontsize = 11.
        clabel_fontsize = 9.5
        axis_tick_fontsize = 9.5
        keep_frames = True
    
        print('\nCreating and saving contour plots ...\n')
        
        get_rounded_str = contourplots.utils.get_rounded_str
        
        # corners of triangular infeasible region
        infeas_ll = 100*yields[0], titers[0]
        infeas_ul = 100*yields[0], titers[-1]
        
        infeas_ur_titer_index = -1
        top_titer_indicator_array = indicator_array[:, -1, :][0]
        infeas_ur_yield_index = np.where(np.isnan(top_titer_indicator_array))[0][-1] + 1
        infeas_ur = 100*yields[infeas_ur_yield_index], titers[infeas_ur_titer_index]
        
        #%% MPSP plot
        
        # MPSP_w_levels, MPSP_w_ticks, MPSP_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
        MPSP_w_levels = np.arange(0., 8.25, 0.25)
        MPSP_cbar_ticks = np.arange(0., 8.1, 1.)
        MPSP_w_ticks = [1., 1.25, 1.5, 1.75, 2., 2.5,  3., 4, 5, 6, 8,]
        # MPSP_w_levels = np.arange(0., 15.5, 0.5)
        
        contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=indicator_array_for_plot, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                        x_data=yields_for_plot, # x axis values
                                        # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                        y_data=titers_for_plot, # y axis values
                                        z_data=productivities, # z axis values
                                        x_label=x_label, # title of the x axis
                                        y_label=y_label, # title of the y axis
                                        z_label=z_label, # title of the z axis
                                        w_label=MPSP_w_label, # title of the color axis
                                        x_ticks=x_ticks,
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
                                        animated_contourplot_filename=filename+'_MPSP_y_t_sims', # file name to save animated contourplot as (no extensions)
                                        keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                        axis_title_fonts=axis_title_fonts,
                                        clabel_fontsize = clabel_fontsize,
                                        default_fontsize = default_fontsize,
                                        axis_tick_fontsize = axis_tick_fontsize,
                                        # comparison_range=market_range,
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
                                        
                                        add_shapes = {
                                            # coords as tuple of tuples: (color, zorder),
                                            (infeas_ll, infeas_ur, infeas_ul): ('white', 2), # infeasible region smoothing
                                            },
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" (",] * 4,
                                        units_closing_brackets = [")",] * 4,
                                        # label_over_color='white',
                                        keep_gifs=False,
                                        )
        
        contourplots.animated_contourplot(
                                        # w_data_vs_x_y_at_multiple_z=[[[MPSP_y_t_f(y,t,ms_fit_coeffs, cs_fit_coeffs) for y in yields] for t in titers]], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                        
                                        w_data_vs_x_y_at_multiple_z=fit_indicator_array_for_plot,
                                        x_data=yields_for_plot, # x axis values
                                        # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                        y_data=titers_for_plot, # y axis values
                                        z_data=productivities, # z axis values
                                        x_label=x_label, # title of the x axis
                                        y_label=y_label, # title of the y axis
                                        z_label=z_label, # title of the z axis
                                        w_label=MPSP_w_label, # title of the color axis
                                        x_ticks=x_ticks,
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
                                        animated_contourplot_filename=filename+'_MPSP_y_t_fit', # file name to save animated contourplot as (no extensions)
                                        keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                        axis_title_fonts=axis_title_fonts,
                                        clabel_fontsize = clabel_fontsize,
                                        default_fontsize = default_fontsize,
                                        axis_tick_fontsize = axis_tick_fontsize,
                                        # comparison_range=market_range,
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
                                        
                                        add_shapes = {
                                            # coords as tuple of tuples: (color, zorder),
                                            (infeas_ll, infeas_ur, infeas_ul): ('white', 2), # infeasible region smoothing
                                            },
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" (",] * 4,
                                        units_closing_brackets = [")",] * 4,
                                        # label_over_color='white',
                                        keep_gifs=False,
                                        )
        
        #%% Relative impact plot
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        # rel_impact_arr = np.array([[10*ti/yi for yi in ys_plot] for ti in ts_plot])
        # rel_impact_arr[np.where(rel_impact_arr>11.)] = 11.
        ax = plt.subplot()
        im = ax.contourf(ys_plot, ts_plot, 
                         rel_impact_arr)
        # ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
        # plt.xlabel('Chlorite dismutase initial conc. [mol/L]')
        # plt.ylabel('Perchlorate reductase initial conc. [mol/L]')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
    
        # ax.clabel(im, fmt='%3.0f', colors='black', fontsize=12)
    
        plt.colorbar(im, cax=cax, label='Relative impact of yield vs titer')
    
        plt.figure()
        
        rel_impact_w_levels = np.arange(0., 10.2, 0.25)
        rel_impact_cbar_ticks = np.arange(0., 10.2, 1.)
        rel_impact_w_ticks = [0, 0.25, 0.5, 1., 2., 4., 8., 10.]
        # rel_impact_w_levels = np.arange(0., 15.5, 0.5)
        
        rel_impact_arr[np.where(np.isnan(indicator_array_for_eval[0]))] = np.nan
        contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=[rel_impact_arr], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                        x_data=ys_plot, # x axis values
                                        # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                        y_data=ts_plot, # y axis values
                                        z_data=productivities, # z axis values
                                        x_label=x_label, # title of the x axis
                                        y_label=y_label, # title of the y axis
                                        z_label=z_label, # title of the z axis
                                        w_label=rel_impact_w_label, # title of the color axis
                                        x_ticks=100*x_ticks,
                                        y_ticks=y_ticks,
                                        z_ticks=z_ticks,
                                        w_levels=rel_impact_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                        w_ticks=rel_impact_w_ticks, # labeled, lined contours; a subset of w_levels
                                        x_units=x_units,
                                        y_units=y_units,
                                        z_units=z_units,
                                        w_units=rel_impact_units,
                                        # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                        fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                        cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                        cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                        extend_cmap='max',
                                        cbar_ticks=rel_impact_cbar_ticks,
                                        z_marker_color='g', # default matplotlib color names
                                        fps=fps, # animation frames (z values traversed) per second
                                        n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                        animated_contourplot_filename=filename+'_rel_impact_y_t_sims', # file name to save animated contourplot as (no extensions)
                                        keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                        axis_title_fonts=axis_title_fonts,
                                        clabel_fontsize = clabel_fontsize,
                                        default_fontsize = default_fontsize,
                                        axis_tick_fontsize = axis_tick_fontsize,
                                        # comparison_range=market_range,
                                        n_minor_ticks = 1,
                                        cbar_n_minor_ticks = 3,
                                        # manual_clabels_regular = {
                                        #     rel_impact_w_ticks[5]: (45,28),
                                        #     },
                                        # additional_points ={(73, 62.5):('D', 'w', 6)},
                                        fill_bottom_with_cmap_over_color=False, # for TRY
                                        # bottom_fill_bounds = ((0,0), 
                                        #                       (5,60.),
                                        #                       (95,10.)),
                                        # # zoom_data_scale=5,
                                        # text_boxes = {'>4.00': [(5,5), 'white']},
                                        
                                        add_shapes = {
                                            # coords as tuple of tuples: (color, zorder),
                                            (infeas_ll, infeas_ur, infeas_ul): ('white', 2), # infeasible region smoothing
                                            },
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" (",] * 4,
                                        units_closing_brackets = [")",] * 4,
                                        # label_over_color='white',
                                        keep_gifs=False,
                                        )
        
    #%% Validate against external data
    if external_data_fit_and_plot:
        folder_name = 'external_data_for_validation'
        product_ID = 'product'
        print('\nValidate with external data:\n')
        plt.rcParams['font.sans-serif'] = "Arial Unicode"
        plt.rcParams['font.size'] = str(12)

        #%% https://doi.org/10.1002/jctb.7690
        # Figure in Jankovic et al. 2024(above DOI) showing data from
        # Jankovic et al. 2023 for ethanol (https://doi.org/10.1016/j.seppur.2023.124320 ) and 
        # Jankovic et al. 2023 for isobutanol (https://doi.org/10.1016/j.ecmx.2023.100520 )
        
        # Tools: Aspen Plus, Aspen Plus
        # Products: ethanol, isobutanol
        # Feedstocks: crude ethanol, crude isobutanol (only modeled separations)
        # Separations: 
            # ethanol: vacuum distillation + extractive distillation (ethylene glycol)
            # isobutanol: Vacuum distillation and a novel hybrid combination of gas stripping and vacuum evaporation were coupled with atmospheric azeotropic distillation
        
        # 6 points, 4 points
        print('\n#1: Janković et al., normalized annual cost of separation ($/kg) vs titer (wt%)')
        filename = folder_name+'/'+'jankovic_et_al_2023_ethanol_7a'
        product_ID = 'product'
        
        ext_data1 = pd.read_csv(filename+'.csv')
        ts_ext1 = ext_data1['titer (wt%)']
        mpsps_ext1 = ext_data1['annual cost of separation ($/kg)']
        
        m,c,r = fit_shifted_rect_hyperbola(ts_ext1, mpsps_ext1)
        
        fig = plt.figure()
        ax = plt.subplot(111)
        
        x_ticks = np.arange(0.0, 0.071, 0.01)
        y_ticks = np.arange(0., 0.401, 0.05)
        
        format_ax(ax, x_ticks, y_ticks,
                  )
        
        # ax = plt.subplot(411)
        
        # ax.set_ylabel(f'MPSP [$/kg {product_ID}] * titer[g/L]')
        # ax.set_xlabel(f'{product_ID} titer [g/L]')
        # ax.scatter(p_curr, [p_curr*arr], label='simulated')
        # ax.plot(p_curr, [m*t + c for t in p_curr], label='fit') # MPSP * titer vs titer straight line fit
        
        # ax = plt.subplot(412)
        ax.set_ylabel(f'MPSP [$/kg {product_ID}]')
        ax.set_xlabel(f'{product_ID} titer [g/L]')
        ax.scatter(ts_ext1, mpsps_ext1, label='simulated')
        ts_ext1_arr = np.array(ts_ext1)
        ts_ext1_fit_plot = np.linspace(ts_ext1_arr[0], ts_ext1_arr[-1], 100)
        ax.plot(ts_ext1_fit_plot, [m + c/t for t in ts_ext1_fit_plot], label='fit') # resulting MPSP vs titer curve
        # plt.legend(loc=(
        #         np.linspace(x_ticks[0], x_ticks[-1], 100)[-4], # x-coord
        #         np.linspace(y_ticks[0], y_ticks[-1], 100)[-20],# y-coord
        #         ))
        plt.legend(edgecolor ='gray')
        
        # ax = plt.subplot(413)
        # ax.set_ylabel(f'Residuals [$/kg {product_ID}]')
        # ax.set_xlabel(f'{product_ID} titer [g/L]')
        # ax.scatter(p_curr, arr - np.array([m + c/t for t in p_curr]))
        
        # ax = plt.subplot(414)
        # ax.set_ylabel('Standardized residuals [-]')
        # ax.set_xlabel(f'{product_ID} titer [g/L]')
        # ax.scatter(p_curr, (arr - np.array([m + c/t for t in p_curr]))/arr)
        # plt.show()
        Rsq = get_Rsq(np.array(mpsps_ext1), 
                      np.array([shifted_rect_hyperbola_f(p, m, c) for p in ts_ext1])
                      )
        print(Rsq)
        print(m, c)
        
        textstr = "$\mathrm{R}^{2}$" + " = " + "%.3f"%(Rsq)
        props = dict(boxstyle='round', 
                     # facecolor='wheat', 
                     facecolor='white', 
                     alpha=0.5,
                     # edgecolor ='gray',
                     )
        
        ax.text(
                np.linspace(x_ticks[0], x_ticks[-1], 100)[-32], # x-coord
                # np.linspace(x_ticks[0], x_ticks[-1], 100)[15], # x-coord
                np.linspace(y_ticks[0], y_ticks[-1], 100)[-30], # y-coord
                textstr, 
                # transform=ax.transAxes, 
                fontsize=12,
                # verticalalignment='top', 
                bbox=props)
        
        plt.savefig(filename+'.png', dpi=300, 
                    bbox_inches='tight',
                    # facecolor=plt.get_facecolor(),
                    # transparent=False,
                    )
        
        # for p_curr, arr in zip(p_currs, arrs):
        #     p_curr = np.array(p_curr)
        #     arr = np.array(arr)
            
        #     # m,c,r = fit_straight_line(p_curr, p_curr*arr)
        #     m,c,r = fit_shifted_rect_hyperbola(p_curr, arr)
            
        #     fig = plt.figure()
        #     ax = plt.subplot(111)
            
        #     # ax = plt.subplot(411)
            
        #     # ax.set_ylabel(f'MPSP [$/kg {product_ID}] * titer[g/L]')
        #     # ax.set_xlabel(f'{product_ID} titer [g/L]')
        #     # ax.scatter(p_curr, [p_curr*arr], label='simulated')
        #     # ax.plot(p_curr, [m*t + c for t in p_curr], label='fit') # MPSP * titer vs titer straight line fit
            
        #     # ax = plt.subplot(412)
        #     ax.set_ylabel(f'MPSP [$/kg {product_ID}]')
        #     ax.set_xlabel(f'{product_ID} titer [g/L]')
        #     ax.scatter(p_curr, arr, label='simulated')
        #     ax.plot(p_curr, [m + c/t for t in p_curr], label='fit') # resulting MPSP vs titer curve
        #     plt.legend()
            
        #     # ax = plt.subplot(413)
        #     # ax.set_ylabel(f'Residuals [$/kg {product_ID}]')
        #     # ax.set_xlabel(f'{product_ID} titer [g/L]')
        #     # ax.scatter(p_curr, arr - np.array([m + c/t for t in p_curr]))
            
        #     # ax = plt.subplot(414)
        #     # ax.set_ylabel('Standardized residuals [-]')
        #     # ax.set_xlabel(f'{product_ID} titer [g/L]')
        #     # ax.scatter(p_curr, (arr - np.array([m + c/t for t in p_curr]))/arr)
        #     # plt.show()
            
        #     print(get_Rsq(np.array(arr), 
        #                   np.array([shifted_rect_hyperbola_f(p, m, c) for p in p_curr])
        #                   )
        #           )
        #     print(m, c)
        #%%
        def titer_f(y, MPSP, a, b, c, d):
            inv_MPSP_minus_a = 1/(MPSP-a)
            y_minus_h = y - b*inv_MPSP_minus_a
            # y_minus_h[np.where(y_minus_h<=0)]=1e-6
            return A(MPSP, a, b, c, d) / y_minus_h\
                    + c*inv_MPSP_minus_a
        
        def A(MPSP, a, b, c, d):
            inv_MPSP_minus_a = 1/(MPSP-a)
            return (d+b*c*inv_MPSP_minus_a)*inv_MPSP_minus_a

        #%% https://doi.org/10.1021/acssuschemeng.7b01729
        # Tool: SuperPro
        # Feedstock: corn
        # Product: 3-HP (anaerobic ferm)
        # Separation: solvent extraction (methyl isobutyl ketone) & distillation
        print('\n#2: Gunukula et al., minimum selling price ($/kg) vs yield (g/g) and titer (g/L)')
        # 31 y-t combinations
        filename = folder_name+'/'+'gunukula_2017_SI_3-HP_b'
        ext_data2 = pd.read_csv(filename+'.csv')
        ys_ext2, ts_ext2 = ext_data2['yield (g/g)'], ext_data2['titer (g/L)']
        mpsps_ext2 = ext_data2['MPSP ($/kg)']
        a_ext2, b_ext2, c_ext2, d_ext2, r = fit_shifted_rect_hyperbola_two_param((ys_ext2, ts_ext2), 
                                                                                   mpsps_ext2,)
        Rsq = get_Rsq(np.array(mpsps_ext2), 
                      np.array([shifted_rect_hyperbola_two_param(y, t, a_ext2, b_ext2, c_ext2, d_ext2)
                                for y,t in zip(ys_ext2, ts_ext2)]))
        print(Rsq)
        fig, axs = plt.subplots(1, 2, constrained_layout=True)
        
        ax = axs[0]
        
        x_ticks = np.arange(0., 1.01, 0.2)
        y_ticks = np.arange(0., 200.01, 20)
        
        format_ax(ax, x_ticks, y_ticks,
                  )
        
        ax.set_xlabel(f'Yield [g/g]')
        ax.set_ylabel(f'{product_ID} Titer [g/L]')
        
        mpsps_ext2 = np.array(mpsps_ext2)
        MPSP_indices = np.where(mpsps_ext2[:-1] != mpsps_ext2[1:])[0]
        MPSP_indices += 1
        # MPSP_indices = [0] + list(MPSP_indices)
        
        curr_MPSP_ind = 0
        next_MPSP_ind = MPSP_indices[0]
        
        plot_this_iter = True
        for i, MPSP in zip(range(len(mpsps_ext2)), mpsps_ext2):
            
            curr_MPSP = mpsps_ext2[curr_MPSP_ind]
            ys_curr = ys_ext2[curr_MPSP_ind:next_MPSP_ind]
            ts_curr = ts_ext2[curr_MPSP_ind:next_MPSP_ind]
            
            # print('\n')
            # print(np.array(ys_curr))
            # print(np.array(ts_curr))
            
            if plot_this_iter:
                ax.plot(ys_curr, 
                        ts_curr,
                        # label='simulated',
                        label=str(round(curr_MPSP, 2)),
                        color='black')
                plot_this_iter = False
                
            if not next_MPSP_ind is None:
                curr_MPSP_ind = 1*next_MPSP_ind
                try:
                    next_MPSP_ind = MPSP_indices[list(MPSP_indices).index(curr_MPSP_ind)+1]
                except:
                    next_MPSP_ind=None
                plot_this_iter = True
        
        labelLines(ax.get_lines(), zorder=2.5)
        
        
        ax = axs[1]
        format_ax(ax, x_ticks, y_ticks,
                  )
        
        ax.set_xlabel(f'Yield [g/g]')
        ax.set_ylabel(f'{product_ID} Titer [g/L]')
        
        curr_MPSP_ind = 0
        next_MPSP_ind = MPSP_indices[0]
        
        plot_this_iter = True
        for i, MPSP in zip(range(len(mpsps_ext2)), mpsps_ext2):
            # ys_curr = ys_ext2[curr_MPSP_ind:next_MPSP_ind]
            # ts_fit = [titer_f(y, MPSP, a_ext2, b_ext2, c_ext2, d_ext2) for y in ys_curr]
            
            
            # print('\n')
            # print(np.array(ys_curr))
            # print(ts_fit)
            
            if plot_this_iter:
                curr_MPSP = mpsps_ext2[curr_MPSP_ind]
                # MPSP = mpsps_ext2[curr_MPSP_ind]
                h = b_ext2/(curr_MPSP-a_ext2)
                k = c_ext2/(curr_MPSP-a_ext2)
                
                
                ys_curr = np.linspace(h*1.001, 1, 200)
                ts_fit = np.array([titer_f(y, curr_MPSP, a_ext2, b_ext2, c_ext2, d_ext2) for y in ys_curr])
                inds_plot = np.where(ts_fit>0)[0]
                # print(ts_fit)
                # print('plot')
                ax.plot(ys_curr[inds_plot], 
                        ts_fit[inds_plot],
                        # label='simulated',
                        label=str(round(curr_MPSP, 2)),
                        color='black')
                plot_this_iter = False
            
            # print('ind', curr_MPSP_ind, next_MPSP_ind)
            if not next_MPSP_ind is None:
                curr_MPSP_ind = 1*next_MPSP_ind
                try:
                    next_MPSP_ind = MPSP_indices[list(MPSP_indices).index(curr_MPSP_ind)+1]
                except:
                    next_MPSP_ind=None
                plot_this_iter = True
        
        labelLines(ax.get_lines(), zorder=2.5)
        
        
        textstr = "$\mathrm{R}^{2}$" + " = " + "%.3f"%(Rsq)
        props = dict(boxstyle='round', 
                     # facecolor='wheat', 
                     facecolor='white', 
                     alpha=1,
                     # edgecolor ='gray',
                     )
        
        ax.text(
                np.linspace(x_ticks[0], x_ticks[-1], 100)[-45], # x-coord
                # np.linspace(x_ticks[0], x_ticks[-1], 100)[15], # x-coord
                np.linspace(y_ticks[0], y_ticks[-1], 100)[-1]*1.05, # y-coord
                textstr, 
                # transform=ax.transAxes, 
                fontsize=12,
                # verticalalignment='top', 
                bbox=props)
        
        plt.savefig(filename+'.png', dpi=300, 
                    bbox_inches='tight',
                    # facecolor=plt.get_facecolor(),
                    # transparent=False,
                    )
            
        # Tool: SuperPro
        # Feedstock: corn
        # Product: 1,3-propanediol (aerobic ferm)
        # Separation: solvent extraction & distillation
        
        
        # Tool: SuperPro
        # Feedstock: corn
        # Product: adipic acid
        # Separation: (microfiltration, centrifugation), distillation, ultrafiltration, crystallization, ion exchange, drying
        
        
        # Tool: SuperPro
        # Feedstock: corn
        # Product: succinic acid
        # Separation: (microfiltration, centrifugation), distillation, ultrafiltration, crystallization, ion exchange, drying
        
        
        # Tool: SuperPro
        # Feedstock: corn
        # Product: 1,3-propanediol (microaerobic ferm)
        # Separation: solvent extraction & distillation
        
        # Tool: SuperPro
        # Feedstock: corn
        # Product: isobutanol
        # Separation: distillation, decantation, distillation
        
        #%% https://doi.org/10.1007/s10098-024-02843-w
        # Tool: Aspen Plus
        # Feedstock: sugarcane A-molasses
        # Product: 2,3-BDO
        # Separation: solvent extraction (oleyl alcohol) & distillation
        print('\n#3: Sikazwe et al., minimum selling price ($/kg) vs yield (g/g) and titer (g/L)')
        # 11 y-t combinations
        filename = folder_name+'/'+'Sikazwe_2023_5_b'
        ext_data2 = pd.read_csv(filename+'.csv')
        ys_ext3, ts_ext3 = ext_data2['yield (g/g)'], ext_data2['titer (g/L)']
        mpsps_ext3 = ext_data2['MPSP ($/kg)']/1000
        a_ext3, b_ext3, c_ext3, d_ext3, r = fit_shifted_rect_hyperbola_two_param((ys_ext3, ts_ext3), 
                                                                                   mpsps_ext3,)
        Rsq = get_Rsq(np.array(mpsps_ext3), 
                      np.array([shifted_rect_hyperbola_two_param(y, t, a_ext3, b_ext3, c_ext3, d_ext3)
                                for y,t in zip(ys_ext3, ts_ext3)]))
        print(Rsq)
        fig, axs = plt.subplots(1, 2, constrained_layout=True)
        
        ax = axs[0]
        
        x_ticks = np.arange(0.3, 0.501, 0.05)
        y_ticks = np.arange(50., 150.01, 25)
        
        format_ax(ax, x_ticks, y_ticks,
                  )
        
        ax.set_xlabel(f'Yield [g/g]')
        ax.set_ylabel(f'{product_ID} Titer [g/L]')
        
        mpsps_ext3 = np.array(mpsps_ext3)
        MPSP_indices = np.where(mpsps_ext3[:-1] != mpsps_ext3[1:])[0]
        MPSP_indices += 1
        # MPSP_indices = [0] + list(MPSP_indices)
        
        curr_MPSP_ind = 0
        next_MPSP_ind = MPSP_indices[0]
        
        plot_this_iter = True
        for i, MPSP in zip(range(len(mpsps_ext3)), mpsps_ext3):
            
            
            curr_MPSP = mpsps_ext2[curr_MPSP_ind]
            ys_curr = ys_ext3[curr_MPSP_ind:next_MPSP_ind]
            ts_curr = ts_ext3[curr_MPSP_ind:next_MPSP_ind]
            
            # print('\n')
            # print(np.array(ys_curr))
            # print(np.array(ts_curr))
            
            if plot_this_iter:
                ax.plot(ys_curr, 
                        ts_curr, 
                        # label='simulated',
                        label = str(round(curr_MPSP, 2)),
                        color='black')
                plot_this_iter = False
            if not next_MPSP_ind is None:
                curr_MPSP_ind = 1*next_MPSP_ind
                try:
                    next_MPSP_ind = MPSP_indices[list(MPSP_indices).index(curr_MPSP_ind)+1]
                except:
                    next_MPSP_ind=None
                plot_this_iter = True
            
        
        labelLines(ax.get_lines(), zorder=2.5)
        
        
        ax = axs[1]
        format_ax(ax, x_ticks, y_ticks,
                  )
        
        ax.set_xlabel(f'Yield [g/g]')
        ax.set_ylabel(f'{product_ID} Titer [g/L]')
        
        curr_MPSP_ind = 0
        next_MPSP_ind = MPSP_indices[0]
        
        plot_this_iter = True
        for i, MPSP in zip(range(len(mpsps_ext3)), mpsps_ext3):
            # ys_curr = ys_ext3[curr_MPSP_ind:next_MPSP_ind]
            # ts_fit = [titer_f(y, MPSP, a_ext3, b_ext3, c_ext3, d_ext3) for y in ys_curr]
            
            
            # print('\n')
            # print(np.array(ys_curr))
            # print(ts_fit)
            
            if plot_this_iter:
                h = b_ext3/(MPSP-a_ext3)
                k = c_ext3/(MPSP-a_ext3)
                MPSP = mpsps_ext3[curr_MPSP_ind]
                
                ys_curr = np.linspace(h*1.001, 1, 200)
                ts_fit = np.array([titer_f(y, MPSP, a_ext3, b_ext3, c_ext3, d_ext3) for y in ys_curr])
                inds_plot = np.where(ts_fit>0)[0]
                # print(ts_fit)
                # print('plot')
                ax.plot(ys_curr[inds_plot], 
                        ts_fit[inds_plot], 
                        # label='simulated',
                        label = str(round(curr_MPSP, 2)),
                        color='black')
                plot_this_iter = False
            
            # print('ind', curr_MPSP_ind, next_MPSP_ind)
            if not next_MPSP_ind is None:
                curr_MPSP_ind = 1*next_MPSP_ind
                try:
                    next_MPSP_ind = MPSP_indices[list(MPSP_indices).index(curr_MPSP_ind)+1]
                except:
                    next_MPSP_ind=None
                plot_this_iter = True
        
        labelLines(ax.get_lines(), zorder=2.5)
        
        
        textstr = "$\mathrm{R}^{2}$" + " = " + "%.3f"%(Rsq)
        props = dict(boxstyle='round', 
                     # facecolor='wheat', 
                     facecolor='white', 
                     alpha=1,
                     # edgecolor ='gray',
                     )
        
        ax.text(
                np.linspace(x_ticks[0], x_ticks[-1], 100)[-45], # x-coord
                # np.linspace(x_ticks[0], x_ticks[-1], 100)[15], # x-coord
                np.linspace(y_ticks[0], y_ticks[-1], 100)[-1]*1.05, # y-coord
                textstr, 
                # transform=ax.transAxes, 
                fontsize=12,
                # verticalalignment='top', 
                bbox=props)
        
        plt.savefig(filename+'.png', dpi=300, 
                    bbox_inches='tight',
                    # facecolor=plt.get_facecolor(),
                    # transparent=False,
                    )
    
    #%% Print coefficients again for convenience
    print('\n')
    print(filename)
    print(a_fit, b_fit, c_fit, d_fit)
    print(a_fit*true_g, b_fit*true_g, c_fit*true_g, d_fit*true_g)
    