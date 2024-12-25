# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 11:42:11 2024

@author: sarangbhagwat
"""

import numpy as np
import pandas as pd
import itertools
import pickle
import os
from fermentation_insights.utils import fit_shifted_rect_hyperbola_two_param, get_feasible_TY_samples

np.random.seed(4153)

#%% Load baseline TRY
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

product = product_ID = 'HP'
# additional_tag = '0.5x_baselineprod'
additional_tag = 'neutral'
# additional_tag = ''
feedstock = 'cornstover'

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

# theoretical max fermentation yield, if loaded yields are in %theoretical max
theo_max_yield = None # kg-fermentation-product-produced / kg-glucose-eq.
if product in ('HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol'): theo_max_yield = 1.
elif product in ('TAL', 'TAL_SA'): theo_max_yield = 0.4667
elif product in ('succinic', 'succinic_neutral'): theo_max_yield = 1.311

yields*=theo_max_yield
inflection_product_yields*=theo_max_yield

true_g = recovery

#%% Import uncertainty analysis stuff

from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
import contourplots
import biosteam as bst
print('\n\nLoading system ...')
from biorefineries import HP

models = None

if feedstock=='glucose':
    from biorefineries.HP.models.glucose import models_glucose_improved_separations
    models = models_glucose_improved_separations
elif feedstock=='sugarcane':
    from biorefineries.HP.models.sugarcane import models_sc_improved_separations
    models = models_sc_improved_separations
elif feedstock=='corn':
    from biorefineries.HP.models.corn import models_corn_improved_separations
    models = models_corn_improved_separations
elif feedstock=='cornstover':
    from biorefineries.HP.models.cornstover import models_cs_improved_separations
    models = models_cs_improved_separations

print('\nLoaded system.')
from datetime import datetime
from biosteam.utils import TicToc
import os

from biorefineries.HP.analyses.full.plot_utils import plot_kde_formatted
from matplotlib.colors import hex2color

chdir = os.chdir
HP_filepath = HP.__file__.replace('\\__init__.py', '')
HP_results_filepath = HP_filepath + '\\analyses\\results\\'
model = models.HP_model

system = HP_sys = models.HP_sys
spec = models.spec
unit_groups = models.unit_groups

tea = models.HP_tea
lca = models.HP_lca
get_adjusted_MSP = models.get_adjusted_MSP
run_bugfix_barrage = models.run_bugfix_barrage

# per_kg_AA_to_per_kg_SA = models.per_kg_AA_to_per_kg_SA

f = bst.main_flowsheet
u, s = f.unit, f.stream

# %% 

N_simulations_per_TRY_combo = 500 # 6000

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]

notification_interval = 10

feedstock_tag = feedstock
product_tag = 'Acrylic' if product=='HP' else product

mode = '300L_FGI'

scenario_name = mode

product_folder = 'acrylic_acid_product' if product_tag=='Acrylic' else 'HP_salt_product'

parameter_distributions_filename = HP_filepath+\
        f'\\analyses\\full\\parameter_distributions\\{product_folder}\\'+\
        f'parameter-distributions_{feedstock_tag}_{product_tag}_' + mode + '.xlsx' 

#%% Bugfix barrage (without set_production_capacity)
baseline_spec = {'spec_1': spec.baseline_yield,
                 'spec_2': spec.baseline_titer,
                 'spec_3': spec.baseline_productivity,}

system=model._system
def reset_and_reload():
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    try:
        print('Resetting cache and emptying recycles ...')
        system.reset_cache()
        system.empty_recycles()
        print('Loading and simulating with baseline specifications ...')
        spec.load_specifications(**baseline_spec)
        # spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
        # spec.set_production_capacity(spec.desired_annual_production)
        system.simulate()
        print('Loading and simulating with required specifications ...')
        spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
        # spec.set_production_capacity(spec.desired_annual_production)
        system.simulate()
    except:
        spec.spec_1, spec.spec_2, spec.spec_3 = spec_1, spec_2, spec_3
        spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
        # spec.set_production_capacity(spec.desired_annual_production)
        system.simulate()
        
    
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    # spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    # spec.set_production_capacity(spec.desired_annual_production)
    system.simulate()

F403 = u.F403
def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        if 'length' in str(e).lower() or 'in subtract' in str(e).lower() or 'in log' in str(e).lower():
            raise(e)
            system.reset_cache()
            system.empty_recycles()
            F403.heat_utilities = []
            F403._V_first_effect = 0.144444
            F403.run()
            F403._design()
            F403.simulate()
            # spec.set_production_capacity(spec.desired_annual_production)
            # spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
            system.simulate()
        else:
            try:
                reset_and_switch_solver('fixedpoint')
            except Exception as e:
                print(str(e))
                try:
                    reset_and_switch_solver('aitken')
                except Exception as e:
                    print(str(e))
                    # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                    print("Bugfix barrage failed.\n")
                    # breakpoint()
                    raise e

#%% Model specification (without set_production_capacity)
pre_fermenter_units_path = list(spec.reactor.get_upstream_units())
pre_fermenter_units_path.reverse()
def model_specification():
    # for i in system.products: i.empty()
    try:
        for i in pre_fermenter_units_path: i.simulate()
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        # spec.set_production_capacity(spec.desired_annual_production)
        model._system.simulate()

    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        # breakpoint()
        # raise e
        if 'sugar concentration' in str_e:
            # flowsheet('AcrylicAcid').F_mass /= 1000.
            raise e
        elif 'length' in str(e).lower() or 'in subtract' in str(e).lower() or 'in log' in str(e).lower():
            raise e
        else:
            # breakpoint()
            run_bugfix_barrage()
            
model.specification = model_specification
spec.reactor.neutralization = False if not 'neutral' in product+additional_tag else True
model.specification()
print(get_adjusted_MSP())

#%% Create TRY combinations 

steps = [5, 5]
# ys = np.linspace(yields[0], yields[-1], steps[0])
# ts = np.linspace(titers[0], titers[-1], steps[1])

def MPSP_sim_f(y, t):
    # spec.load_specifications(y, t, spec.spec_3)
    spec.spec_1, spec.spec_2 = y, t
    spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    model._system.simulate()
    return get_adjusted_MSP()

yts = get_feasible_TY_samples(yields, titers, steps, MPSP_sim_f, theo_max_yield=theo_max_yield)
yts = np.array(yts)

spec.load_specifications(**baseline_spec)
model.specification()

#%%
# overall_results_dict = {(y,t): {'Baseline':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}, 
#                             # 'GWP Breakdown':{}, 'FEC Breakdown':{},
#                             },
#                 'Uncertainty':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
#                 'Sensitivity':{'Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
#                                'p-val Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}}},}
#                 for (y,t) in list(itertools.product(ys,ts))}

overall_results_dict = {}

#%% Load parameter distributions and samples
print(f'\n\nLoading parameter distributions ({mode}) ,,,')
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename, models.namespace_dict)

# load_additional_params()
print(f'\nLoaded parameter distributions ({mode}).')

parameters = model.get_parameters()

print('\n\nLoading samples ...')
samples = model.sample(N=N_simulations_per_TRY_combo, rule='L')
model.load_samples(samples)
print('\nLoaded samples.')

# ## Change working directory to biorefineries\\HP\\analyses\\results
# chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##

model.exception_hook = 'warn'


#%% Run uncertainty analyses across TRY

errors_dict = {}

for yt in yts:
    try:
        y, t = yt
        print('\n\n------------------------------------------------------------------------------------')
        print(f'\nPerforming uncertainty analysis at yield = {np.round(y,2)} and titer = {np.round(t,2)} ...')
        spec.spec_1 = y/theo_max_yield
        spec.spec_2 = t
        # try:
            # results_dict = overall_results_dict[(y,t)]
        overall_results_dict[(y,t)] = results_dict = {'Baseline':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}, 
                                    # 'GWP Breakdown':{}, 'FEC Breakdown':{},
                                    },
                                    'Uncertainty':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}, 'Recovery':{}},
                                    'Sensitivity':{'Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}},
                                                    'p-val Spearman':{'MPSP':{}, 'GWP100a':{}, 'FEC':{}}},}
        
        
        model.load_samples(samples)
        # Initial baseline simulation
        print('\n\nSimulating baseline ...')
        
        baseline_initial = model.metrics_at_baseline()
        # spec.set_production_capacity()
        baseline_initial = model.metrics_at_baseline()
        
        # baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
        #                         columns=baseline_initial.keys())
        
        results_dict['Baseline']['MPSP'][mode] = get_adjusted_MSP()
        results_dict['Baseline']['GWP100a'][mode] = tot_GWP = lca.GWP
        results_dict['Baseline']['FEC'][mode] = tot_FEC = lca.FEC
    
        print(f"\nSimulated baseline. MPSP = ${round(results_dict['Baseline']['MPSP'][mode],2)}/kg.")
        
        print('\n\nEvaluating ...')
        try:
            model.evaluate(notify=notification_interval, autoload=None, autosave=None, file=None)
        except Exception as e:
            print('\n'+str(e))
            breakpoint()
            # raise(e)
            
        print('\nFinished evaluation.')
        
        
        # Final baseline simulation
        print('\n\nRe-simulating baseline ...')
        
        baseline_end = model.metrics_at_baseline()
        # spec.set_production_capacity()
        baseline_end = model.metrics_at_baseline()
        
        print(f"\nRe-simulated baseline. MPSP = ${round(get_adjusted_MSP(),2)}/kg.")
        
        
        table = model.table
        model.table = model.table.dropna()
    
        # results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
        results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price [$/kg AA]']
        results_dict['Uncertainty']['GWP100a'][mode] = model.table.Biorefinery['Total gwp100a [kg-CO2-eq/kg]'] # GWP or gwp
        results_dict['Uncertainty']['FEC'][mode] = model.table.Biorefinery['Total FEC [MJ/kg]']
        results_dict['Uncertainty']['Recovery'][mode] = model.table.Biorefinery['Product recovery [%]']
        df_rho, df_p = model.spearman_r()
        
        # results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
        results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg AA]']
        results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = df_rho['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
        results_dict['Sensitivity']['Spearman']['FEC'][mode] = df_rho['Biorefinery', 'Total FEC [MJ/kg]']
        
        # results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
        results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price [$/kg AA]']
        results_dict['Sensitivity']['p-val Spearman']['GWP100a'][mode] = df_p['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
        results_dict['Sensitivity']['p-val Spearman']['FEC'][mode] = df_p['Biorefinery', 'Total FEC [MJ/kg]']
        
        
        # Parameters
        parameters = model.get_parameters()
        index_parameters = len(model.get_baseline_sample())
        parameter_values = model.table.iloc[:, :index_parameters].copy()
        
        results_dict['Parameters'] = parameter_values
    except:
        print(f'\n\nFailed evaluation at y = {y}, t = {t}.')
    # except Exception as e:
    #     print(str(e))
    #     errors_dict[(y,t)] = e

#%% Save results for uncertainty analyses across TRY

with open(filename+'_ord_uncertainty.pkl', 'wb') as f:
    pickle.dump(overall_results_dict, f)

with open(filename+'_ord_uncertainty.pkl', 'rb') as f:
    overall_results_dict = pickle.load(f)
    
#%% Utils for MPSP_f

from fermentation_insights.utils import shifted_rect_hyperbola_two_param

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


#%% Estimate coeffs at each set of 25 TRY combos for N_simulations_per_TRY_combo alternative sets of parameter values

ap, bp, cp, dp, g = [], [], [], [], []

final_product_MW = 72.06266

fermentation_product_MW = 90.07794

for i in range(N_simulations_per_TRY_combo):
    MPSPs = []
    recoveries = []
    yts_fit = []
    for yt in yts:
        y, t = yt
        try:
            results_dict = overall_results_dict[(y,t)]
            recoveries.append(results_dict['Uncertainty']['Recovery'][mode][i] * final_product_MW/fermentation_product_MW)
            MPSPs.append(results_dict['Uncertainty']['MPSP'][mode][i])
            yts_fit.append((y,t))
        except: # error, likely infeasible region'
            print(f'Error for {yt}.')
            # pass
        
    try:
        yts_fit = np.array(yts_fit)
        a_fit, b_fit, c_fit, d_fit, r = fit_shifted_rect_hyperbola_two_param((yts_fit[:,0], yts_fit[:,1]), 
                                                                            MPSPs,
                                                                            [1,1,1,1],
                                                                            # [0,0,0,0,0],
                                                                            )
        print(get_Rsq(np.array(MPSPs), 
                      np.array([shifted_rect_hyperbola_two_param(y, t, a_fit, b_fit, c_fit, d_fit) 
                                  for y,t in yts_fit])))
        ap.append(a_fit)
        bp.append(b_fit)
        cp.append(c_fit)
        dp.append(d_fit)
        g.append(np.median(np.array(recoveries)))
    except Exception as e:
        print(f'Error for {i}: {str(e)}')
        ap.append(np.nan)
        bp.append(np.nan)
        cp.append(np.nan)
        dp.append(np.nan)
        g.append(np.nan)
    

a = np.array(ap)*g
b = np.array(bp)*g
c = np.array(cp)*g
d = np.array(dp)*g
g = np.array(g)

#%%
print([np.median(i[np.where(~np.isnan(i))]) for i in [a, b, c, d, g]])
print([np.percentile(i[np.where(~np.isnan(i))], 5) for i in [a, b, c, d, g]])
print([np.percentile(i[np.where(~np.isnan(i))], 95) for i in [a, b, c, d, g]])

#%% Spearman's -- coeffs w.r.t. uncertain parameters
df_dict = {'a':a,
 'b':b,
 'c':c,
 'd':d,
 'g':g}

df_dict.update({k:v for k,v in overall_results_dict[(y,t)]['Parameters'].items()})
df = pd.DataFrame.from_dict(df_dict)
spearman = df.corr(method='spearman')

#%% Save values of parameters and (fit) coefficients

with open(filename+'_dfd_coeffs_uncertainty.pkl', 'wb') as f:
    pickle.dump(df_dict, f)

#%% 

with open(filename+'_dfd_coeffs_uncertainty.pkl', 'rb') as f:
    df_dict = pickle.load(f)
    
#%%
# spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
# model.metrics_at_baseline()

# MPSPs = []
# recoveries = []
# yts = []

# for y in ys:
#     for t in ts:
#         spec.spec_1 = y
#         spec.spec_2 = t
#         model.specification()
#         MPSPs.append(model.metrics[3]())
#         recoveries.append(model.metrics[5]()* final_product_MW/fermentation_product_MW)
#         yts.append((y,t))
        
# yts = np.array(yts)
# a_fit, b_fit, c_fit, d_fit, r = fit_shifted_rect_hyperbola_two_param((yts[:,0], yts[:,1]), 
#                                                                     MPSPs,
#                                                                     [1,1,1,1],
#                                                                     # [0,0,0,0,0],
#                                                                     )
# g_fit = np.median(np.array(recoveries))

# print(a_fit*g_fit, b_fit*g_fit, c_fit*g_fit, d_fit*g_fit)