# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 23:03:53 2025

@author: sarangbhagwat
"""
#%% Initial
import numpy as np
import os
import multiprocessing
import itertools

keep_varnames = [
                 'feedstocks',
                 'neutralizations',
                 'productivity_lavels',
                 'steps',
                 'os',
                 'np'
                 'run_TRY_analysis', 
                 'productivities', 
                 'inflection_product_yields', 
                 'product_ID',]

feedstocks = ['glucose', 'sugarcane', 'corn', 'cornstover']
neutralizations = [False, True]
productivity_labels = [None, '0.2bp', '5.0bp']
steps = 3

processes = 12

def get_productivities(baseline_productivity,
                       multipliers=np.array([1., 0.2, 5.])):
    return multipliers*baseline_productivity

#%% succinic
product_ID = 'succinic'
inflection_product_yields = np.array([1-0.186])
productivities = get_productivities(baseline_productivity=0.6573)
from biorefineries.succinic.analyses.TRY_analysis_FGI import run_TRY_analysis

for feedstock in feedstocks:
    for neutralization in neutralizations:
        for prod, prod_label in zip(productivities, productivity_labels):
            print(f'\n\n{product_ID}, {neutralization}, {prod_label}, {feedstock}')
            yields, titers,\
            results_metric_1, results_metric_2, results_metric_3,\
            results_metric_4 = run_TRY_analysis(feedstock=feedstock, neutralization=neutralization,
                                                return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
                                                productivities=np.array([prod]))
    
            
            tag = product_ID + '_'
            if neutralization:
                tag += 'neutral_'
            if prod_label is not None:
                tag += prod_label + '_'
            tag += feedstock
            
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
            np.save(f'{tag}_productivities', productivities)
            
            del yields, titers, inflection_product_yields
            del results_metric_1, results_metric_2, results_metric_3, results_metric_4
            
            import gc
            gc.collect()

#%% TAL
product_ID = 'TAL'
productivities = get_productivities(baseline_productivity=0.12)
from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_FGI import run_TRY_analysis

all_combs = itertools.product(feedstocks, productivities)


def func(comb):
    feedstock, prod = comb
    return run_TRY_analysis(feedstock=feedstock,
            return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
            productivities=np.array([prod]),
            which_fig='insights',
            print_status_every_n_simulations=500,)

print(f'\n\n{product_ID}\n\n')

results = None
with multiprocessing.Pool(processes) as pool:
    results = pool.map(func, all_combs)

for r, prod_label in zip(results, productivity_labels):
    yields, titers, inflection_product_yields,\
    results_metric_1, results_metric_2, results_metric_3,\
    results_metric_4, results_metric_5, results_metric_6 = r
    tag = product_ID + '_'
    if prod_label is not None:
        tag += prod_label + '_'
    tag += feedstock
    
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
    np.save(f'{tag}_productivities', productivities)
  
    del yields, titers, inflection_product_yields
    del results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6

    import gc
    gc.collect()

del results

import gc
gc.collect()
        
#%% TAL_SA
product_ID = 'TAL_SA'
productivities = get_productivities(baseline_productivity=0.12)
from biorefineries.TAL.analyses.fermentation.TRY_analysis_SA_FGI import run_TRY_analysis

for feedstock in feedstocks:
    # for neutralization in neutralizations:
    for prod, prod_label in zip(productivities, productivity_labels):
        print(f'\n\n{product_ID}, {prod_label}, {feedstock}')
        yields, titers, inflection_product_yields,\
        results_metric_1, results_metric_2, results_metric_3,\
        results_metric_4, results_metric_5, results_metric_6 = run_TRY_analysis(feedstock=feedstock,
                                                               return_YT_arrays=True, plot=False, steps_for_default_YT=steps,
                                                               productivities=np.array([prod]),
                                                               which_fig='insights',
                                                               print_status_every_n_simulations=500,)
        
        tag = product_ID + '_'
        if prod_label is not None:
            tag += prod_label + '_'
        tag += feedstock
        
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
        np.save(f'{tag}_productivities', productivities)
      
        del yields, titers, inflection_product_yields
        del results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6
        
        import gc
        gc.collect()

