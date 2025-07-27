# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 11:42:11 2024

@author: sarangbhagwat
"""

from numba import config
# config.DISABLE_JIT = True
import numpy as np
import pandas as pd
import itertools
import pickle
import os
from fermentation_insights.utils import fit_shifted_rect_hyperbola_two_param, get_feasible_TY_samples
from biosteam.units import MultiEffectEvaporator

def run_uncertainties_coeffs(product, additional_tag, feedstock, steps_TRY=5, N_simulations_per_TRY_combo=1000, notification_interval=100):
    import numpy as np
    np.random.seed(4153)
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    
    product_ID = product
    # additional_tag = '0.5x_baselineprod'
    # additional_tag = ''
    feedstock = feedstock
    
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
    import biosteam as bst
    print('\n\nLoading system ...')
    
    models = None
    model = None
    system = None
    spec = None
    tea, lca = None, None
    get_adjusted_MSP = None
    unit_groups = None
    
    biorefinery_filepath = None
    biorefinery_results_filepath = None
    
    if product=='TAL':
        from biorefineries import TAL
        
        if feedstock=='glucose':
            from biorefineries.TAL.models import models_TAL_solubility_exploit_glucose
            models = models_TAL_solubility_exploit_glucose
        elif feedstock=='sugarcane':
            from biorefineries.TAL.models import models_TAL_solubility_exploit_sugarcane
            models = models_TAL_solubility_exploit_sugarcane
        elif feedstock=='corn':
            from biorefineries.TAL.models import models_TAL_solubility_exploit_corn
            models = models_TAL_solubility_exploit_corn
        elif feedstock=='cornstover':
            from biorefineries.TAL.models import models_TAL_solubility_exploit_cornstover
            models = models_TAL_solubility_exploit_cornstover
        
        model = models.TAL_model
        system = TAL_sys = models.TAL_sys
        spec = models.spec
        tea = models.TAL_tea
        lca = models.TAL_lca
        get_adjusted_MSP = models.get_adjusted_MSP
        unit_groups = models.unit_groups
        
        biorefinery_filepath = TAL.__file__.replace('\\__init__.py', '')
        biorefinery_results_filepath = biorefinery_filepath + '\\analyses\\results\\'
    
    elif product=='TAL_SA':
        from biorefineries import TAL
        
        if feedstock=='glucose':
            from biorefineries.TAL.models import models_SA_solubility_exploit_glucose
            models = models_SA_solubility_exploit_glucose
        elif feedstock=='sugarcane':
            from biorefineries.TAL.models import models_SA_solubility_exploit_sugarcane
            models = models_SA_solubility_exploit_sugarcane
        elif feedstock=='corn':
            from biorefineries.TAL.models import models_SA_solubility_exploit_corn
            models = models_SA_solubility_exploit_corn
        elif feedstock=='cornstover':
            from biorefineries.TAL.models import models_SA_solubility_exploit_cornstover
            models = models_SA_solubility_exploit_cornstover
        
        model = models.TAL_model
        system = TAL_sys = models.TAL_sys
        spec = models.spec
        unit_groups = models.unit_groups
        tea = models.TAL_tea
        lca = models.TAL_lca
        get_adjusted_MSP = models.get_adjusted_MSP
        run_bugfix_barrage = models.run_bugfix_barrage
        
        biorefinery_filepath = TAL.__file__.replace('\\__init__.py', '')
        biorefinery_results_filepath = biorefinery_filepath + '\\analyses\\results\\'
        
    elif product=='HP' and additional_tag=='':
        from biorefineries import HP
        
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
        
        model = models.HP_model
        system = HP_sys = models.HP_sys
        spec = models.spec
        unit_groups = models.unit_groups
        tea = models.HP_tea
        lca = models.HP_lca
        get_adjusted_MSP = models.get_adjusted_MSP
        run_bugfix_barrage = models.run_bugfix_barrage
        
        biorefinery_filepath = HP.__file__.replace('\\__init__.py', '')
        biorefinery_results_filepath = biorefinery_filepath + '\\analyses\\results\\'
        
    elif product=='HP' and ('hexanol' in additional_tag):
        from biorefineries import HP
        
        if feedstock=='glucose':
            from biorefineries.HP.models.glucose import models_glucose_hexanol
            models = models_glucose_hexanol
        elif feedstock=='sugarcane':
            from biorefineries.HP.models.sugarcane import models_sc_hexanol
            models = models_sc_hexanol
        elif feedstock=='corn':
            from biorefineries.HP.models.corn import models_corn_hexanol
            models = models_corn_hexanol
        elif feedstock=='cornstover':
            from biorefineries.HP.models.cornstover import models_cs_hexanol
            models = models_cs_hexanol
        
        model = models.HP_model
        system = HP_sys = models.HP_sys
        spec = models.spec
        unit_groups = models.unit_groups
        tea = models.HP_tea
        lca = models.HP_lca
        get_adjusted_MSP = models.get_adjusted_MSP
        biorefinery_filepath = HP.__file__.replace('\\__init__.py', '')
        biorefinery_results_filepath = biorefinery_filepath + '\\analyses\\results\\'
        
    elif product=='succinic':
        from biorefineries import succinic
        
        if feedstock=='glucose':
            from biorefineries.succinic import models_glucose
            models = models_glucose
        elif feedstock=='sugarcane':
            from biorefineries.succinic import models_corn
            models = models_corn
        elif feedstock=='corn':
            from biorefineries.succinic import models_sugarcane
            models = models_sugarcane
        elif feedstock=='cornstover':
            from biorefineries.succinic import models_cornstover
            models = models_cornstover
        
        model = models.succinic_model
        system = succinic_sys = models.succinic_sys
        spec = models.spec
        unit_groups = models.unit_groups
        tea = models.succinic_tea
        lca = models.succinic_LCA
        get_adjusted_MSP = models.get_adjusted_MSP
        run_bugfix_barrage = models.run_bugfix_barrage
        
        biorefinery_filepath = succinic.__file__.replace('\\__init__.py', '')
        biorefinery_results_filepath = biorefinery_filepath + '\\analyses\\results\\'
        
    print('\nLoaded system.')
    print(get_adjusted_MSP())
    chdir = os.chdir
    
    run_bugfix_barrage = models.run_bugfix_barrage
    
    # per_kg_AA_to_per_kg_SA = models.per_kg_AA_to_per_kg_SA
    
    f = bst.main_flowsheet
    u, s = f.unit, f.stream
    
    # %% 
    N_simulations_per_TRY_combo = N_simulations_per_TRY_combo # 6000
    
    percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]
    
    feedstock_tag = feedstock
    product_tag = product.replace('_neutral', '')
    
    print('\nFetching parameter distributions ...')
    
    if product=='TAL_SA':
        product_tag = 'Sorbate'
    elif product in ('HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol'):
        product_tag = 'Acrylic'
    
    mode = None
    parameter_distributions_filename = None
    
    if product in ('TAL', 'TAL_SA'):
        mode = 'A_FGI'
        parameter_distributions_filename = biorefinery_filepath+\
                f'\\analyses\\full\\parameter_distributions\\'+\
                f'parameter-distributions_{product_tag}_' + mode + f'_{feedstock_tag}'+'.xlsx' 
        
    elif product in ('HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol'):
        mode = '300L_FGI'
        product_folder = 'acrylic_acid_product' if product_tag=='Acrylic' else 'HP_salt_product'

        parameter_distributions_filename = biorefinery_filepath+\
                f'\\analyses\\full\\parameter_distributions\\{product_folder}\\'+\
                f'parameter-distributions_{feedstock_tag}_{product_tag}_' + mode + '.xlsx' 
    
    elif product in ('succinic', 'succinic_neutral'):
        mode = f'pilot-scale_batch_{feedstock_tag}_FGI' 
        parameter_distributions_filename = biorefinery_filepath+'\\analyses\\parameter_distributions\\'+'parameter-distributions_' + mode + '.xlsx'

    print('\nFetched parameter distributions.')
    
    #%% Bugfix barrage (without set_production_capacity)
    baseline_spec = {'spec_1': spec.baseline_yield,
                     'spec_2': spec.baseline_titer,
                     'spec_3': spec.baseline_productivity,}
    
    system=model._system
    def reset_V_first_effect_evaporators(system):
        for unit in system.units:
            if isinstance(unit, MultiEffectEvaporator): 
                unit._V_first_effect = None
    
    def run_evaporators(system):
        for unit in system.units:
            if isinstance(unit, MultiEffectEvaporator): 
                try:
                    unit._run()
                except:
                    pass
                finally:
                    unit._run()
    
    def reset_and_reload():
        spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
        try:
            print('Resetting cache and emptying recycles ...')
            system.reset_cache()
            system.empty_recycles()
            run_evaporators(system)
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
        run_evaporators(system)
        system.converge_method = solver_ID
        print(f"Trying {solver_ID} ...")
        # spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        # spec.set_production_capacity(spec.desired_annual_production)
        system.simulate()
    
    # F403 = u.F403
    def traditional_bugfix_sequence():
        try:
            reset_and_reload()
        except Exception as e:
            print(str(e))
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
                    system.simulate()
                    # breakpoint()
                    raise e
                finally:
                    system.simulate()
                    
    def run_bugfix_barrage(e):
        if 'multieffectevaporator' in str(e).lower():
            try:
                print('Resetting _V_first_effect of all evaporators to None ...')
                reset_V_first_effect_evaporators(system)
                print('Running all evaporators ...')
                run_evaporators(system)
                print('Simulating ...')
                system.simulate()
                
            except:
                traditional_bugfix_sequence()
        else:
            traditional_bugfix_sequence()
    ###############################
    
    #%% Model specification (without set_production_capacity)
    pre_fermenter_units_path = list(spec.reactor.get_upstream_units())
    pre_fermenter_units_path.reverse()
    def model_specification():
        try:
            for i in pre_fermenter_units_path: i.simulate()
            spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
            # spec.set_production_capacity(spec.desired_annual_production)
            # system.simulate()
            model._system.simulate()
        
        except Exception as e:
            str_e = str(e).lower()
            print('Error in model spec: %s'%str_e)
            # raise e
            if 'sugar concentration' in str_e:
                # flowsheet('AcrylicAcid').F_mass /= 1000.
                raise e
            elif 'length' in str(e).lower() or 'in subtract' in str(e).lower() or 'in log' in str(e).lower():
                raise e
            else:
                if product=='succinic':
                    try:
                        model._system.simulate()
                        for i in pre_fermenter_units_path: i.simulate()
                        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
                        model._system.simulate()
                    except:
                        run_bugfix_barrage(e)
                else:
                    run_bugfix_barrage(e)
                    
    model.specification = model_specification
    
    #%% Create TRY combinations 
    
    steps = [steps_TRY, steps_TRY]
    # ys = np.linspace(yields[0], yields[-1], steps[0])
    # ts = np.linspace(titers[0], titers[-1], steps[1])
    
    def MPSP_sim_f(y, t):
        # spec.load_specifications(y, t, spec.spec_3)
        spec.spec_1, spec.spec_2 = y, t
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        model._system.simulate()
        return get_adjusted_MSP()
    
    print('\nSimulating to get feasible TY samples ...')
    yts = get_feasible_TY_samples(yields, titers, steps, MPSP_sim_f, theo_max_yield=theo_max_yield)
    yts = np.array(yts)
    print(f'\nFeasible TY samples obtained: {yts}')
    
    print('\nSimulating baseline ...')
    spec.load_specifications(**baseline_spec)
    model.specification()
    print('\nSimulated baseline.')
    
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
    
    model.exception_hook = 'warn'
    
    #%% Run uncertainty analyses across TRY
    
    errors_dict = {}
    
    n_unc_analyses = len(yts)
    for i, yt in zip(range(n_unc_analyses), yts):
        
        y, t = yt
        print('\n\n------------------------------------------------------------------------------------')
        print(f'\nPerforming uncertainty analysis ({i}/{n_unc_analyses}), at yield = {np.round(y,2)} and titer = {np.round(t,2)} ...')
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
            # print('\n'+str(e))
            # breakpoint()
            raise(e)
            
        print('\nFinished evaluation.')
        
        
        # Final baseline simulation
        print('\n\nRe-simulating baseline ...')
        
        baseline_end = model.metrics_at_baseline()
        # spec.set_production_capacity()
        baseline_end = model.metrics_at_baseline()
        
        print(f"\nRe-simulated baseline. MPSP = ${round(get_adjusted_MSP(),2)}/kg.")
        
        
        table = model.table
        model.table = model.table.dropna()
    
        if product=='TAL':
            # results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price [$/kg TAL]']
            results_dict['Uncertainty']['GWP100a'][mode] = model.table.Biorefinery['Total gwp100a [kg-CO2-eq/kg]'] # GWP or gwp
            results_dict['Uncertainty']['FEC'][mode] = model.table.Biorefinery['Total FEC [MJ/kg]']
            results_dict['Uncertainty']['Recovery'][mode] = model.table.Biorefinery['Product recovery [%]']
            df_rho, df_p = model.spearman_r()
            
            # results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg TAL]']
            results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = df_rho['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
            results_dict['Sensitivity']['Spearman']['FEC'][mode] = df_rho['Biorefinery', 'Total FEC [MJ/kg]']
            
            # results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price [$/kg TAL]']
            results_dict['Sensitivity']['p-val Spearman']['GWP100a'][mode] = df_p['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
            results_dict['Sensitivity']['p-val Spearman']['FEC'][mode] = df_p['Biorefinery', 'Total FEC [MJ/kg]']
        
        elif product=='TAL_SA':
            # results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price [$/kg KSA]']
            results_dict['Uncertainty']['GWP100a'][mode] = model.table.Biorefinery['Total gwp100a [kg-CO2-eq/kg]'] # GWP or gwp
            results_dict['Uncertainty']['FEC'][mode] = model.table.Biorefinery['Total FEC [MJ/kg]']
            results_dict['Uncertainty']['Recovery'][mode] = model.table.Biorefinery['Product recovery [%]']
            df_rho, df_p = model.spearman_r()
            
            # results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg KSA]']
            results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = df_rho['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
            results_dict['Sensitivity']['Spearman']['FEC'][mode] = df_rho['Biorefinery', 'Total FEC [MJ/kg]']
            
            # results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price [$/kg KSA]']
            results_dict['Sensitivity']['p-val Spearman']['GWP100a'][mode] = df_p['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
            results_dict['Sensitivity']['p-val Spearman']['FEC'][mode] = df_p['Biorefinery', 'Total FEC [MJ/kg]']
            
        elif product=='HP' and additional_tag=='':
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
        
        elif product=='HP' and ('hexanol' in additional_tag):
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
            
        elif product=='succinic':
            # results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Uncertainty']['MPSP'][mode] = model.table.Biorefinery['Adjusted minimum selling price [$/kg]']
            results_dict['Uncertainty']['GWP100a'][mode] = model.table.Biorefinery['Total gwp100a [kg-CO2-eq/kg]'] # GWP or gwp
            results_dict['Uncertainty']['FEC'][mode] = model.table.Biorefinery['Total FEC [MJ/kg]']
            results_dict['Uncertainty']['Recovery'][mode] = model.table.Biorefinery['Product recovery [%]']
            df_rho, df_p = model.spearman_r()
            
            # results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Sensitivity']['Spearman']['MPSP'][mode] = df_rho['Biorefinery', 'Adjusted minimum selling price [$/kg]']
            results_dict['Sensitivity']['Spearman']['GWP100a'][mode] = df_rho['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
            results_dict['Sensitivity']['Spearman']['FEC'][mode] = df_rho['Biorefinery', 'Total FEC [MJ/kg]']
            
            # results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price - as sorbic acid [$/kg SA-eq.]']
            results_dict['Sensitivity']['p-val Spearman']['MPSP'][mode] = df_p['Biorefinery', 'Adjusted minimum selling price [$/kg]']
            results_dict['Sensitivity']['p-val Spearman']['GWP100a'][mode] = df_p['Biorefinery', 'Total gwp100a [kg-CO2-eq/kg]']
            results_dict['Sensitivity']['p-val Spearman']['FEC'][mode] = df_p['Biorefinery', 'Total FEC [MJ/kg]']
            
        # Parameters
        parameters = model.get_parameters()
        index_parameters = len(model.get_baseline_sample())
        parameter_values = model.table.iloc[:, :index_parameters].copy()
        
        results_dict['Parameters'] = parameter_values
            
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
    
    final_product_MW, fermentation_product_MW = None, None
    
    if product=='TAL':
        final_product_MW = 126.11004
        fermentation_product_MW = 126.11004
    elif product=='TAL_SA':
        final_product_MW = 150.21688
        fermentation_product_MW = 126.11004
    elif product=='HP':
        final_product_MW = 72.06266
        fermentation_product_MW = 90.07794
        
    elif product=='succinic':
        final_product_MW = 118.09
        fermentation_product_MW = 118.09
    
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
            g.append(np.median(np.array(np.nan)))
        
    
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