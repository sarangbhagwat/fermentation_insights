# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 23:03:53 2025

@author: sarangbhagwat
"""

#%% Initial
# from numba import config
# config.DISABLE_JIT = True
import numpy as np
import os
import multiprocessing
import itertools

def get_productivities(baseline_productivity,
                       multipliers=np.array([1., 
                                             0.2, 
                                             5.,
                                             ])):
    return multipliers*baseline_productivity

products = [
            # 'TAL',
            # 'TAL_SA', 
            # 'HP', 
            # 'HP_hexanol', 
            'succinic',
            ]
feedstocks = ['glucose', 
              # 'sugarcane', 
              # 'corn', 
              # 'cornstover'
              ]

neutralizations_all = {
    'TAL': [False],
    'TAL_SA': [False],
    'HP': [False, True],
    'HP_hexanol': [False, True],
    'succinic': [False, True],
    }

productivity_labels = [None, '0.2bp', '5.0bp']

productivities_all = {
    'TAL': get_productivities(baseline_productivity=0.12),
    'TAL_SA': get_productivities(baseline_productivity=0.12),
    'HP': get_productivities(baseline_productivity=0.276),
    'HP_hexanol': get_productivities(baseline_productivity=0.276),
    'succinic': get_productivities(baseline_productivity=0.6573)}

steps = 50

processes = 32

from biorefineries.succinic.analyses.TRY_analysis_FGI import run_TRY_analysis as run_TRY_analysis_succinic
from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_FGI import run_TRY_analysis as run_TRY_analysis_TAL
from biorefineries.TAL.analyses.fermentation.TRY_analysis_SA_FGI import run_TRY_analysis as run_TRY_analysis_SA
from biorefineries.HP.analyses.fermentation.TRY_analysis_FGI import run_TRY_analysis as run_TRY_analysis_HP
from biorefineries.HP.analyses.fermentation.TRY_analysis_hexanol_FGI import run_TRY_analysis as run_TRY_analysis_HP_hexanol

all_combs = []

for product in products:
    all_combs += list(itertools.product([product], neutralizations_all[product], feedstocks, productivities_all[product]))

def func(comb):
    product, neutralization, feedstock, prod = comb
    
    if product=='TAL':
        return run_TRY_analysis_TAL(feedstock=feedstock,
                return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
                productivities=np.array([prod]),
                which_fig='insights',
                print_status_every_n_simulations=500,)
    elif product=='TAL_SA':
        return run_TRY_analysis_SA(feedstock=feedstock,
                return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
                productivities=np.array([prod]),
                which_fig='insights',
                print_status_every_n_simulations=500,)
    elif product=='succinic':
        return run_TRY_analysis_succinic(feedstock=feedstock, neutralization=neutralization,
                return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
                productivities=np.array([prod]))
    elif product=='HP':
        return run_TRY_analysis_HP(feedstock=feedstock, neutralization=neutralization,
                return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
                productivities=np.array([prod]))
    elif product=='HP_hexanol':
        return run_TRY_analysis_HP_hexanol(feedstock=feedstock, neutralization=neutralization,
                return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
                productivities=np.array([prod]))
    
# print(f'\n\n{product}\n\n')

def save_results(comb, r):
    product, neutralization, feedstock, prod = comb
    prod_label = productivity_labels[list(productivities_all[product]).index(prod)]
    
    tag = product + '_'
    if neutralization:
        tag += 'neutral_'
    if prod_label is not None:
        tag += prod_label + '_'
    tag += feedstock
    
    tag = tag.replace('HP_hexanol_neutral', 'HP_neutral_hexanol')
    
    if product=='TAL' or product=='TAL_SA':
        yields, titers, inflection_product_yields,\
        results_metric_1, results_metric_2, results_metric_3,\
        results_metric_4, results_metric_5, results_metric_6 = r
    
        os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
        
        np.save(f'{tag}_MPSP', results_metric_1)
        np.save(f'{tag}_GWP', results_metric_2)
        np.save(f'{tag}_FEC', results_metric_3)
        np.save(f'{tag}_AOC', results_metric_5)
        np.save(f'{tag}_TCI', results_metric_6)
        np.save(f'{tag}_recoveries', results_metric_4)
        np.save(f'{tag}_inflection_product_yields', inflection_product_yields)
        np.save(f'{tag}_yields', yields)
        np.save(f'{tag}_titers', titers)
        np.save(f'{tag}_productivities', productivities_all[product])
        
        print(f'\nSaved results for {tag}')
    
    elif product=='succinic':
        yields, titers, inflection_product_yields,\
        results_metric_1, results_metric_2, results_metric_3,\
        results_metric_4 = r
    
        os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
        
        np.save(f'{tag}_MPSP', results_metric_1)
        np.save(f'{tag}_GWP', results_metric_2)
        np.save(f'{tag}_FEC', results_metric_3)
        # np.save('{tag}_AOC', results_metric_4)
        # np.save('{tag}_TCI', results_metric_5)
        np.save(f'{tag}_recoveries', results_metric_4)
        np.save(f'{tag}_inflection_product_yields', inflection_product_yields)
        np.save(f'{tag}_yields', yields)
        np.save(f'{tag}_titers', titers)
        np.save(f'{tag}_productivities', productivities_all[product])
        
        print(f'\nSaved results for {tag}')
        
    elif product=='HP' or product=='HP_hexanol':
        yields, titers, inflection_product_yields,\
        results_metric_1, results_metric_2, results_metric_3,\
        results_metric_4, results_metric_5, results_metric_6 = r
    
        os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
        
        np.save(f'{tag}_MPSP', results_metric_1)
        np.save(f'{tag}_GWP', results_metric_2)
        np.save(f'{tag}_FEC', results_metric_3)
        np.save(f'{tag}_AOC', results_metric_4)
        np.save(f'{tag}_TCI', results_metric_5)
        np.save(f'{tag}_recoveries', results_metric_6)
        np.save(f'{tag}_inflection_product_yields', inflection_product_yields)
        np.save(f'{tag}_yields', yields)
        np.save(f'{tag}_titers', titers)
        np.save(f'{tag}_productivities', productivities_all[product])
        
        print(f'\nSaved results for {tag}')
        
def func_with_save(comb):
    r = func(comb)
    save_results(comb, r)
        
#%%
run = False
if __name__ == '__main__' and run:  
    with multiprocessing.Pool(processes) as pool:
        results = pool.map(func, all_combs)
    for r, comb in zip(results, all_combs):
        product, neutralization, feedstock, prod = comb
        prod_label = productivity_labels[list(productivities_all[product]).index(prod)]
        save_results(comb, r)
