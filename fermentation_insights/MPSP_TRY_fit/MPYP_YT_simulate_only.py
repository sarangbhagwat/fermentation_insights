# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 17:09:07 2024

@author: sarangbhagwat
"""

#%% TAL

#%% TAL_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'TAL'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_glucose import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec,\
                                                                TAL_maximum_viable_market_range as market_range
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_glucose_MPSP', results_metric_1)
    np.save('TAL_glucose_GWP', results_metric_2)
    np.save('TAL_glucose_FEC', results_metric_3)
    np.save('TAL_glucose_AOC', results_metric_5)
    np.save('TAL_glucose_TCI', results_metric_6)
    np.save('TAL_glucose_recoveries', results_metric_4)
    np.save('TAL_glucose_inflection_product_yields', inflection_product_yields)
    np.save('TAL_glucose_yields', yields)
    np.save('TAL_glucose_titers', titers)
    np.save('TAL_glucose_productivities', productivities)
else:
    results_metric_1 = np.load('TAL_glucose_MPSP')
    results_metric_2 = np.load('TAL_glucose_GWP')
    results_metric_3 = np.load('TAL_glucose_FEC')
    results_metric_5 = np.load('TAL_glucose_AOC')
    results_metric_6 = np.load('TAL_glucose_TCI')
    inflection_product_yields = np.load('TAL_glucose_inflection_product_yields')

#%% TAL_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'TAL'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_sugarcane import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec,\
                                                                TAL_maximum_viable_market_range as market_range
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_sugarcane_MPSP', results_metric_1)
    np.save('TAL_sugarcane_GWP', results_metric_2)
    np.save('TAL_sugarcane_FEC', results_metric_3)
    np.save('TAL_sugarcane_AOC', results_metric_5)
    np.save('TAL_sugarcane_TCI', results_metric_6)
    np.save('TAL_sugarcane_recoveries', results_metric_4)
    np.save('TAL_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('TAL_sugarcane_yields', yields)
    np.save('TAL_sugarcane_titers', titers)
    np.save('TAL_sugarcane_productivities', productivities)
else:
    results_metric_1 = np.load('TAL_sugarcane_MPSP')
    results_metric_2 = np.load('TAL_sugarcane_GWP')
    results_metric_3 = np.load('TAL_sugarcane_FEC')
    results_metric_5 = np.load('TAL_sugarcane_AOC')
    results_metric_6 = np.load('TAL_sugarcane_TCI')
    inflection_product_yields = np.load('TAL_sugarcane_inflection_product_yields')
            
#%% TAL_corn
run_simulations = True
import numpy as np
import os
product_ID = 'TAL'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_corn import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec,\
                                                                TAL_maximum_viable_market_range as market_range
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_corn_MPSP', results_metric_1)
    np.save('TAL_corn_GWP', results_metric_2)
    np.save('TAL_corn_FEC', results_metric_3)
    np.save('TAL_corn_AOC', results_metric_5)
    np.save('TAL_corn_TCI', results_metric_6)
    np.save('TAL_corn_recoveries', results_metric_4)
    np.save('TAL_corn_inflection_product_yields', inflection_product_yields)
    np.save('TAL_corn_yields', yields)
    np.save('TAL_corn_titers', titers)
    np.save('TAL_corn_productivities', productivities)
else:
    results_metric_1 = np.load('TAL_corn_MPSP')
    results_metric_2 = np.load('TAL_corn_GWP')
    results_metric_3 = np.load('TAL_corn_FEC')
    results_metric_5 = np.load('TAL_corn_AOC')
    results_metric_6 = np.load('TAL_corn_TCI')
    inflection_product_yields = np.load('TAL_corn_inflection_product_yields')
            
#%% TAL_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'TAL'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_TAL_cornstover import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec,\
                                                                TAL_maximum_viable_market_range as market_range
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_cornstover_MPSP', results_metric_1)
    np.save('TAL_cornstover_GWP', results_metric_2)
    np.save('TAL_cornstover_FEC', results_metric_3)
    np.save('TAL_cornstover_AOC', results_metric_5)
    np.save('TAL_cornstover_TCI', results_metric_6)
    np.save('TAL_cornstover_recoveries', results_metric_4)
    np.save('TAL_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('TAL_cornstover_yields', yields)
    np.save('TAL_cornstover_titers', titers)
    np.save('TAL_cornstover_productivities', productivities)
else:
    results_metric_1 = np.load('TAL_cornstover_MPSP')
    results_metric_2 = np.load('TAL_cornstover_GWP')
    results_metric_3 = np.load('TAL_cornstover_FEC')
    results_metric_5 = np.load('TAL_cornstover_AOC')
    results_metric_6 = np.load('TAL_cornstover_TCI')
    inflection_product_yields = np.load('TAL_cornstover_inflection_product_yields')
            

#%% TAL_SA

#%% TAL_SA_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'TAL-SA'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_SA_glucose import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec
                                                                
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_SA_glucose_MPSP', results_metric_1)
    np.save('TAL_SA_glucose_GWP', results_metric_2)
    np.save('TAL_SA_glucose_FEC', results_metric_3)
    np.save('TAL_SA_glucose_AOC', results_metric_5)
    np.save('TAL_SA_glucose_TCI', results_metric_6)
    np.save('TAL_SA_glucose_recoveries', results_metric_4)
    np.save('TAL_SA_glucose_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_glucose_yields', yields)
    np.save('TAL_SA_glucose_titers', titers)
    np.save('TAL_SA_glucose_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_glucose_MPSP')
    results_metric_2 = np.load('TAL_SA_glucose_GWP')
    results_metric_3 = np.load('TAL_SA_glucose_FEC')
    results_metric_5 = np.load('TAL_SA_glucose_AOC')
    results_metric_6 = np.load('TAL_SA_glucose_TCI')
    inflection_product_yields = np.load('TAL_SA_glucose_inflection_product_yields')
            
#%% TAL_SA_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'TAL-SA'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_SA_sugarcane import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec
                                                                
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_SA_sugarcane_MPSP', results_metric_1)
    np.save('TAL_SA_sugarcane_GWP', results_metric_2)
    np.save('TAL_SA_sugarcane_FEC', results_metric_3)
    np.save('TAL_SA_sugarcane_AOC', results_metric_5)
    np.save('TAL_SA_sugarcane_TCI', results_metric_6)
    np.save('TAL_SA_sugarcane_recoveries', results_metric_4)
    np.save('TAL_SA_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_sugarcane_yields', yields)
    np.save('TAL_SA_sugarcane_titers', titers)
    np.save('TAL_SA_sugarcane_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_sugarcane_MPSP')
    results_metric_2 = np.load('TAL_SA_sugarcane_GWP')
    results_metric_3 = np.load('TAL_SA_sugarcane_FEC')
    results_metric_5 = np.load('TAL_SA_sugarcane_AOC')
    results_metric_6 = np.load('TAL_SA_sugarcane_TCI')
    inflection_product_yields = np.load('TAL_SA_sugarcane_inflection_product_yields')
            
#%% TAL_SA_corn
run_simulations = True
import numpy as np
import os
product_ID = 'TAL-SA'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_SA_corn import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec
                                                                
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_SA_corn_MPSP', results_metric_1)
    np.save('TAL_SA_corn_GWP', results_metric_2)
    np.save('TAL_SA_corn_FEC', results_metric_3)
    np.save('TAL_SA_corn_AOC', results_metric_5)
    np.save('TAL_SA_corn_TCI', results_metric_6)
    np.save('TAL_SA_corn_recoveries', results_metric_4)
    np.save('TAL_SA_corn_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_corn_yields', yields)
    np.save('TAL_SA_corn_titers', titers)
    np.save('TAL_SA_corn_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_corn_MPSP')
    results_metric_2 = np.load('TAL_SA_corn_GWP')
    results_metric_3 = np.load('TAL_SA_corn_FEC')
    results_metric_5 = np.load('TAL_SA_corn_AOC')
    results_metric_6 = np.load('TAL_SA_corn_TCI')
    inflection_product_yields = np.load('TAL_SA_corn_inflection_product_yields')
            
#%% TAL_SA_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'TAL-SA'
if run_simulations:
    from biorefineries.TAL.analyses.fermentation.TRY_analysis_SA_cornstover import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, titers_mol_per_mol_total, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                R302, spec
                                                                
    
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_microbe_conversion-R302.regular_citric_acid_conversion,
                                          1-R302.regular_microbe_conversion])
    np.save('TAL_SA_cornstover_MPSP', results_metric_1)
    np.save('TAL_SA_cornstover_GWP', results_metric_2)
    np.save('TAL_SA_cornstover_FEC', results_metric_3)
    np.save('TAL_SA_cornstover_AOC', results_metric_5)
    np.save('TAL_SA_cornstover_TCI', results_metric_6)
    np.save('TAL_SA_cornstover_recoveries', results_metric_4)
    np.save('TAL_SA_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_cornstover_yields', yields)
    np.save('TAL_SA_cornstover_titers', titers)
    np.save('TAL_SA_cornstover_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_cornstover_MPSP')
    results_metric_2 = np.load('TAL_SA_cornstover_GWP')
    results_metric_3 = np.load('TAL_SA_cornstover_FEC')
    results_metric_5 = np.load('TAL_SA_cornstover_AOC')
    results_metric_6 = np.load('TAL_SA_cornstover_TCI')
    inflection_product_yields = np.load('TAL_SA_cornstover_inflection_product_yields')
            
#%% HP

#%% HP_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_glucose_MPSP', results_metric_1)
    np.save('HP_glucose_GWP', results_metric_2)
    np.save('HP_glucose_FEC', results_metric_3)
    np.save('HP_glucose_AOC', results_metric_4)
    np.save('HP_glucose_TCI', results_metric_5)
    np.save('HP_glucose_recoveries', results_metric_6)
    np.save('HP_glucose_inflection_product_yields', inflection_product_yields)
    np.save('HP_glucose_yields', yields)
    np.save('HP_glucose_titers', titers)
    np.save('HP_glucose_productivities', productivities)
else:
    results_metric_1 = np.load('HP_glucose_MPSP')
    results_metric_2 = np.load('HP_glucose_GWP')
    results_metric_3 = np.load('HP_glucose_FEC')
    results_metric_4 = np.load('HP_glucose_AOC')
    results_metric_5 = np.load('HP_glucose_TCI')
    inflection_product_yields = np.load('HP_glucose_inflection_product_yields')

#%% HP_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_sugarcane_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_sugarcane_MPSP', results_metric_1)
    np.save('HP_sugarcane_GWP', results_metric_2)
    np.save('HP_sugarcane_FEC', results_metric_3)
    np.save('HP_sugarcane_AOC', results_metric_4)
    np.save('HP_sugarcane_TCI', results_metric_5)
    np.save('HP_sugarcane_recoveries', results_metric_6)
    np.save('HP_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('HP_sugarcane_yields', yields)
    np.save('HP_sugarcane_titers', titers)
    np.save('HP_sugarcane_productivities', productivities)
else:
    results_metric_1 = np.load('HP_sugarcane_MPSP')
    results_metric_2 = np.load('HP_sugarcane_GWP')
    results_metric_3 = np.load('HP_sugarcane_FEC')
    results_metric_4 = np.load('HP_sugarcane_AOC')
    results_metric_5 = np.load('HP_sugarcane_TCI')
    inflection_product_yields = np.load('HP_sugarcane_inflection_product_yields')

#%% HP_corn
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_corn_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_corn_MPSP', results_metric_1)
    np.save('HP_corn_GWP', results_metric_2)
    np.save('HP_corn_FEC', results_metric_3)
    np.save('HP_corn_AOC', results_metric_4)
    np.save('HP_corn_TCI', results_metric_5)
    np.save('HP_corn_recoveries', results_metric_6)
    np.save('HP_corn_inflection_product_yields', inflection_product_yields)
    np.save('HP_corn_yields', yields)
    np.save('HP_corn_titers', titers)
    np.save('HP_corn_productivities', productivities)
else:
    results_metric_1 = np.load('HP_corn_MPSP')
    results_metric_2 = np.load('HP_corn_GWP')
    results_metric_3 = np.load('HP_corn_FEC')
    results_metric_4 = np.load('HP_corn_AOC')
    results_metric_5 = np.load('HP_corn_TCI')
    inflection_product_yields = np.load('HP_corn_inflection_product_yields')

#%% HP_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_cornstover_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_cornstover_MPSP', results_metric_1)
    np.save('HP_cornstover_GWP', results_metric_2)
    np.save('HP_cornstover_FEC', results_metric_3)
    np.save('HP_cornstover_AOC', results_metric_4)
    np.save('HP_cornstover_TCI', results_metric_5)
    np.save('HP_cornstover_recoveries', results_metric_6)
    np.save('HP_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('HP_cornstover_yields', yields)
    np.save('HP_cornstover_titers', titers)
    np.save('HP_cornstover_productivities', productivities)
else:
    results_metric_1 = np.load('HP_cornstover_MPSP')
    results_metric_2 = np.load('HP_cornstover_GWP')
    results_metric_3 = np.load('HP_cornstover_FEC')
    results_metric_4 = np.load('HP_cornstover_AOC')
    results_metric_5 = np.load('HP_cornstover_TCI')
    inflection_product_yields = np.load('HP_cornstover_inflection_product_yields')

#%% HP_neutral

#%% HP_neutral_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_glucose_MPSP', results_metric_1)
    np.save('HP_neutral_glucose_GWP', results_metric_2)
    np.save('HP_neutral_glucose_FEC', results_metric_3)
    np.save('HP_neutral_glucose_AOC', results_metric_4)
    np.save('HP_neutral_glucose_TCI', results_metric_5)
    np.save('HP_neutral_glucose_recoveries', results_metric_6)
    np.save('HP_neutral_glucose_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_glucose_yields', yields)
    np.save('HP_neutral_glucose_titers', titers)
    np.save('HP_neutral_glucose_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_glucose_MPSP')
    results_metric_2 = np.load('HP_neutral_glucose_GWP')
    results_metric_3 = np.load('HP_neutral_glucose_FEC')
    results_metric_4 = np.load('HP_neutral_glucose_AOC')
    results_metric_5 = np.load('HP_neutral_glucose_TCI')
    inflection_product_yields = np.load('HP_neutral_glucose_inflection_product_yields')

#%% HP_neutral_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_sugarcane_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_sugarcane_MPSP', results_metric_1)
    np.save('HP_neutral_sugarcane_GWP', results_metric_2)
    np.save('HP_neutral_sugarcane_FEC', results_metric_3)
    np.save('HP_neutral_sugarcane_AOC', results_metric_4)
    np.save('HP_neutral_sugarcane_TCI', results_metric_5)
    np.save('HP_neutral_sugarcane_recoveries', results_metric_6)
    np.save('HP_neutral_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_sugarcane_yields', yields)
    np.save('HP_neutral_sugarcane_titers', titers)
    np.save('HP_neutral_sugarcane_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_sugarcane_MPSP')
    results_metric_2 = np.load('HP_neutral_sugarcane_GWP')
    results_metric_3 = np.load('HP_neutral_sugarcane_FEC')
    results_metric_4 = np.load('HP_neutral_sugarcane_AOC')
    results_metric_5 = np.load('HP_neutral_sugarcane_TCI')
    inflection_product_yields = np.load('HP_neutral_sugarcane_inflection_product_yields')

#%% HP_neutral_corn
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_corn_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_corn_MPSP', results_metric_1)
    np.save('HP_neutral_corn_GWP', results_metric_2)
    np.save('HP_neutral_corn_FEC', results_metric_3)
    np.save('HP_neutral_corn_AOC', results_metric_4)
    np.save('HP_neutral_corn_TCI', results_metric_5)
    np.save('HP_neutral_corn_recoveries', results_metric_6)
    np.save('HP_neutral_corn_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_corn_yields', yields)
    np.save('HP_neutral_corn_titers', titers)
    np.save('HP_neutral_corn_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_corn_MPSP')
    results_metric_2 = np.load('HP_neutral_corn_GWP')
    results_metric_3 = np.load('HP_neutral_corn_FEC')
    results_metric_4 = np.load('HP_neutral_corn_AOC')
    results_metric_5 = np.load('HP_neutral_corn_TCI')
    inflection_product_yields = np.load('HP_neutral_corn_inflection_product_yields')

#%% HP_neutral_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_cornstover_Acrylic import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_cornstover_MPSP', results_metric_1)
    np.save('HP_neutral_cornstover_GWP', results_metric_2)
    np.save('HP_neutral_cornstover_FEC', results_metric_3)
    np.save('HP_neutral_cornstover_AOC', results_metric_4)
    np.save('HP_neutral_cornstover_TCI', results_metric_5)
    np.save('HP_neutral_cornstover_recoveries', results_metric_6)
    np.save('HP_neutral_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_cornstover_yields', yields)
    np.save('HP_neutral_cornstover_titers', titers)
    np.save('HP_neutral_cornstover_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_cornstover_MPSP')
    results_metric_2 = np.load('HP_neutral_cornstover_GWP')
    results_metric_3 = np.load('HP_neutral_cornstover_FEC')
    results_metric_4 = np.load('HP_neutral_cornstover_AOC')
    results_metric_5 = np.load('HP_neutral_cornstover_TCI')
    inflection_product_yields = np.load('HP_neutral_cornstover_inflection_product_yields')


#%% HP_hexanol

#%% HP_hexanol_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_hexanol_glucose_MPSP', results_metric_1)
    np.save('HP_hexanol_glucose_GWP', results_metric_2)
    np.save('HP_hexanol_glucose_FEC', results_metric_3)
    np.save('HP_hexanol_glucose_AOC', results_metric_4)
    np.save('HP_hexanol_glucose_TCI', results_metric_5)
    np.save('HP_hexanol_glucose_recoveries', results_metric_6)
    np.save('HP_hexanol_glucose_inflection_product_yields', inflection_product_yields)
    np.save('HP_hexanol_glucose_yields', yields)
    np.save('HP_hexanol_glucose_titers', titers)
    np.save('HP_hexanol_glucose_productivities', productivities)
else:
    results_metric_1 = np.load('HP_hexanol_glucose_MPSP')
    results_metric_2 = np.load('HP_hexanol_glucose_GWP')
    results_metric_3 = np.load('HP_hexanol_glucose_FEC')
    results_metric_4 = np.load('HP_hexanol_glucose_AOC')
    results_metric_5 = np.load('HP_hexanol_glucose_TCI')
    inflection_product_yields = np.load('HP_hexanol_glucose_inflection_product_yields')

#%% HP_hexanol_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_sugarcane_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_hexanol_sugarcane_MPSP', results_metric_1)
    np.save('HP_hexanol_sugarcane_GWP', results_metric_2)
    np.save('HP_hexanol_sugarcane_FEC', results_metric_3)
    np.save('HP_hexanol_sugarcane_AOC', results_metric_4)
    np.save('HP_hexanol_sugarcane_TCI', results_metric_5)
    np.save('HP_hexanol_sugarcane_recoveries', results_metric_6)
    np.save('HP_hexanol_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('HP_hexanol_sugarcane_yields', yields)
    np.save('HP_hexanol_sugarcane_titers', titers)
    np.save('HP_hexanol_sugarcane_productivities', productivities)
else:
    results_metric_1 = np.load('HP_hexanol_sugarcane_MPSP')
    results_metric_2 = np.load('HP_hexanol_sugarcane_GWP')
    results_metric_3 = np.load('HP_hexanol_sugarcane_FEC')
    results_metric_4 = np.load('HP_hexanol_sugarcane_AOC')
    results_metric_5 = np.load('HP_hexanol_sugarcane_TCI')
    inflection_product_yields = np.load('HP_hexanol_sugarcane_inflection_product_yields')

#%% HP_hexanol_corn
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_corn_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_hexanol_corn_MPSP', results_metric_1)
    np.save('HP_hexanol_corn_GWP', results_metric_2)
    np.save('HP_hexanol_corn_FEC', results_metric_3)
    np.save('HP_hexanol_corn_AOC', results_metric_4)
    np.save('HP_hexanol_corn_TCI', results_metric_5)
    np.save('HP_hexanol_corn_recoveries', results_metric_6)
    np.save('HP_hexanol_corn_inflection_product_yields', inflection_product_yields)
    np.save('HP_hexanol_corn_yields', yields)
    np.save('HP_hexanol_corn_titers', titers)
    np.save('HP_hexanol_corn_productivities', productivities)
else:
    results_metric_1 = np.load('HP_hexanol_corn_MPSP')
    results_metric_2 = np.load('HP_hexanol_corn_GWP')
    results_metric_3 = np.load('HP_hexanol_corn_FEC')
    results_metric_4 = np.load('HP_hexanol_corn_AOC')
    results_metric_5 = np.load('HP_hexanol_corn_TCI')
    inflection_product_yields = np.load('HP_hexanol_corn_inflection_product_yields')

#%% HP_hexanol_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_cornstover_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_hexanol_cornstover_MPSP', results_metric_1)
    np.save('HP_hexanol_cornstover_GWP', results_metric_2)
    np.save('HP_hexanol_cornstover_FEC', results_metric_3)
    np.save('HP_hexanol_cornstover_AOC', results_metric_4)
    np.save('HP_hexanol_cornstover_TCI', results_metric_5)
    np.save('HP_hexanol_cornstover_recoveries', results_metric_6)
    np.save('HP_hexanol_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('HP_hexanol_cornstover_yields', yields)
    np.save('HP_hexanol_cornstover_titers', titers)
    np.save('HP_hexanol_cornstover_productivities', productivities)
else:
    results_metric_1 = np.load('HP_hexanol_cornstover_MPSP')
    results_metric_2 = np.load('HP_hexanol_cornstover_GWP')
    results_metric_3 = np.load('HP_hexanol_cornstover_FEC')
    results_metric_4 = np.load('HP_hexanol_cornstover_AOC')
    results_metric_5 = np.load('HP_hexanol_cornstover_TCI')
    inflection_product_yields = np.load('HP_hexanol_cornstover_inflection_product_yields')


#%% HP_hexanol_neutral

#%% HP_neutral_hexanol_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_glucose_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_hexanol_glucose_MPSP', results_metric_1)
    np.save('HP_neutral_hexanol_glucose_GWP', results_metric_2)
    np.save('HP_neutral_hexanol_glucose_FEC', results_metric_3)
    np.save('HP_neutral_hexanol_glucose_AOC', results_metric_4)
    np.save('HP_neutral_hexanol_glucose_TCI', results_metric_5)
    np.save('HP_neutral_hexanol_glucose_recoveries', results_metric_6)
    np.save('HP_neutral_hexanol_glucose_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_hexanol_glucose_yields', yields)
    np.save('HP_neutral_hexanol_glucose_titers', titers)
    np.save('HP_neutral_hexanol_glucose_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_hexanol_glucose_MPSP')
    results_metric_2 = np.load('HP_neutral_hexanol_glucose_GWP')
    results_metric_3 = np.load('HP_neutral_hexanol_glucose_FEC')
    results_metric_4 = np.load('HP_neutral_hexanol_glucose_AOC')
    results_metric_5 = np.load('HP_neutral_hexanol_glucose_TCI')
    inflection_product_yields = np.load('HP_neutral_hexanol_glucose_inflection_product_yields')

#%% HP_neutral_hexanol_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_sugarcane_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_hexanol_sugarcane_MPSP', results_metric_1)
    np.save('HP_neutral_hexanol_sugarcane_GWP', results_metric_2)
    np.save('HP_neutral_hexanol_sugarcane_FEC', results_metric_3)
    np.save('HP_neutral_hexanol_sugarcane_AOC', results_metric_4)
    np.save('HP_neutral_hexanol_sugarcane_TCI', results_metric_5)
    np.save('HP_neutral_hexanol_sugarcane_recoveries', results_metric_6)
    np.save('HP_neutral_hexanol_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_hexanol_sugarcane_yields', yields)
    np.save('HP_neutral_hexanol_sugarcane_titers', titers)
    np.save('HP_neutral_hexanol_sugarcane_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_hexanol_sugarcane_MPSP')
    results_metric_2 = np.load('HP_neutral_hexanol_sugarcane_GWP')
    results_metric_3 = np.load('HP_neutral_hexanol_sugarcane_FEC')
    results_metric_4 = np.load('HP_neutral_hexanol_sugarcane_AOC')
    results_metric_5 = np.load('HP_neutral_hexanol_sugarcane_TCI')
    inflection_product_yields = np.load('HP_neutral_hexanol_sugarcane_inflection_product_yields')

#%% HP_neutral_hexanol_corn
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_corn_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_hexanol_corn_MPSP', results_metric_1)
    np.save('HP_neutral_hexanol_corn_GWP', results_metric_2)
    np.save('HP_neutral_hexanol_corn_FEC', results_metric_3)
    np.save('HP_neutral_hexanol_corn_AOC', results_metric_4)
    np.save('HP_neutral_hexanol_corn_TCI', results_metric_5)
    np.save('HP_neutral_hexanol_corn_recoveries', results_metric_6)
    np.save('HP_neutral_hexanol_corn_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_hexanol_corn_yields', yields)
    np.save('HP_neutral_hexanol_corn_titers', titers)
    np.save('HP_neutral_hexanol_corn_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_hexanol_corn_MPSP')
    results_metric_2 = np.load('HP_neutral_hexanol_corn_GWP')
    results_metric_3 = np.load('HP_neutral_hexanol_corn_FEC')
    results_metric_4 = np.load('HP_neutral_hexanol_corn_AOC')
    results_metric_5 = np.load('HP_neutral_hexanol_corn_TCI')
    inflection_product_yields = np.load('HP_neutral_hexanol_corn_inflection_product_yields')

#%% HP_neutral_hexanol_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'HP'
if run_simulations:
    from biorefineries.HP.analyses.fermentation.TRY_analysis_cornstover_Acrylic_hexanol import results_metric_1, results_metric_2, results_metric_3, results_metric_4, results_metric_5, results_metric_6,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap, get_rounded_str,\
                                                                AA_market_range as market_range,\
                                                                R302, spec, get_AA_MPSP
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    inflection_product_yields = np.array([1-R302.regular_biomass_conversion])
    
    np.save('HP_neutral_hexanol_cornstover_MPSP', results_metric_1)
    np.save('HP_neutral_hexanol_cornstover_GWP', results_metric_2)
    np.save('HP_neutral_hexanol_cornstover_FEC', results_metric_3)
    np.save('HP_neutral_hexanol_cornstover_AOC', results_metric_4)
    np.save('HP_neutral_hexanol_cornstover_TCI', results_metric_5)
    np.save('HP_neutral_hexanol_cornstover_recoveries', results_metric_6)
    np.save('HP_neutral_hexanol_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('HP_neutral_hexanol_cornstover_yields', yields)
    np.save('HP_neutral_hexanol_cornstover_titers', titers)
    np.save('HP_neutral_hexanol_cornstover_productivities', productivities)
else:
    results_metric_1 = np.load('HP_neutral_hexanol_cornstover_MPSP')
    results_metric_2 = np.load('HP_neutral_hexanol_cornstover_GWP')
    results_metric_3 = np.load('HP_neutral_hexanol_cornstover_FEC')
    results_metric_4 = np.load('HP_neutral_hexanol_cornstover_AOC')
    results_metric_5 = np.load('HP_neutral_hexanol_cornstover_TCI')
    inflection_product_yields = np.load('HP_neutral_hexanol_cornstover_inflection_product_yields')


#%% succinic

#%% succinic_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_glucose import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_glucose_MPSP', results_metric_1)
    np.save('succinic_glucose_GWP', results_metric_2)
    np.save('succinic_glucose_FEC', results_metric_3)
    # np.save('succinic_glucose_AOC', results_metric_4)
    # np.save('succinic_glucose_TCI', results_metric_5)
    np.save('succinic_glucose_recoveries', results_metric_4)
    np.save('succinic_glucose_inflection_product_yields', inflection_product_yields)
    np.save('succinic_glucose_yields', yields)
    np.save('succinic_glucose_titers', titers)
    np.save('succinic_glucose_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_glucose_MPSP')
    results_metric_2 = np.load('succinic_glucose_GWP')
    results_metric_3 = np.load('succinic_glucose_FEC')
    recoveries = np.load('succinic_glucose_recoveries')
    # results_metric_4 = np.load('succinic_glucose_AOC')
    # results_metric_5 = np.load('succinic_glucose_TCI')
    inflection_product_yields = np.load('succinic_glucose_inflection_product_yields')

#%% succinic_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_sugarcane import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_sugarcane_MPSP', results_metric_1)
    np.save('succinic_sugarcane_GWP', results_metric_2)
    np.save('succinic_sugarcane_FEC', results_metric_3)
    # np.save('succinic_sugarcane_AOC', results_metric_4)
    # np.save('succinic_sugarcane_TCI', results_metric_5)
    np.save('succinic_sugarcane_recoveries', results_metric_4)
    np.save('succinic_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('succinic_sugarcane_yields', yields)
    np.save('succinic_sugarcane_titers', titers)
    np.save('succinic_sugarcane_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_sugarcane_MPSP')
    results_metric_2 = np.load('succinic_sugarcane_GWP')
    results_metric_3 = np.load('succinic_sugarcane_FEC')
    recoveries = np.load('succinic_sugarcane_recoveries')
    # results_metric_4 = np.load('succinic_sugarcane_AOC')
    # results_metric_5 = np.load('succinic_sugarcane_TCI')
    inflection_product_yields = np.load('succinic_sugarcane_inflection_product_yields')

#%% succinic_corn
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_corn import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_corn_MPSP', results_metric_1)
    np.save('succinic_corn_GWP', results_metric_2)
    np.save('succinic_corn_FEC', results_metric_3)
    # np.save('succinic_corn_AOC', results_metric_4)
    # np.save('succinic_corn_TCI', results_metric_5)
    np.save('succinic_corn_recoveries', results_metric_4)
    np.save('succinic_corn_inflection_product_yields', inflection_product_yields)
    np.save('succinic_corn_yields', yields)
    np.save('succinic_corn_titers', titers)
    np.save('succinic_corn_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_corn_MPSP')
    results_metric_2 = np.load('succinic_corn_GWP')
    results_metric_3 = np.load('succinic_corn_FEC')
    recoveries = np.load('succinic_corn_recoveries')
    # results_metric_4 = np.load('succinic_corn_AOC')
    # results_metric_5 = np.load('succinic_corn_TCI')
    inflection_product_yields = np.load('succinic_corn_inflection_product_yields')

#%% succinic_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_cornstover import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_cornstover_MPSP', results_metric_1)
    np.save('succinic_cornstover_GWP', results_metric_2)
    np.save('succinic_cornstover_FEC', results_metric_3)
    # np.save('succinic_cornstover_AOC', results_metric_4)
    # np.save('succinic_cornstover_TCI', results_metric_5)
    np.save('succinic_cornstover_recoveries', results_metric_4)
    np.save('succinic_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('succinic_cornstover_yields', yields)
    np.save('succinic_cornstover_titers', titers)
    np.save('succinic_cornstover_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_cornstover_MPSP')
    results_metric_2 = np.load('succinic_cornstover_GWP')
    results_metric_3 = np.load('succinic_cornstover_FEC')
    recoveries = np.load('succinic_cornstover_recoveries')
    # results_metric_4 = np.load('succinic_cornstover_AOC')
    # results_metric_5 = np.load('succinic_cornstover_TCI')
    inflection_product_yields = np.load('succinic_cornstover_inflection_product_yields')


#%% succinic_neutral

#%% succinic_neutral_glucose
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_glucose import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_neutral_glucose_MPSP', results_metric_1)
    np.save('succinic_neutral_glucose_GWP', results_metric_2)
    np.save('succinic_neutral_glucose_FEC', results_metric_3)
    # np.save('succinic_neutral_glucose_AOC', results_metric_4)
    # np.save('succinic_neutral_glucose_TCI', results_metric_5)
    np.save('succinic_neutral_glucose_recoveries', results_metric_4)
    np.save('succinic_neutral_glucose_inflection_product_yields', inflection_product_yields)
    np.save('succinic_neutral_glucose_yields', yields)
    np.save('succinic_neutral_glucose_titers', titers)
    np.save('succinic_neutral_glucose_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_neutral_glucose_MPSP')
    results_metric_2 = np.load('succinic_neutral_glucose_GWP')
    results_metric_3 = np.load('succinic_neutral_glucose_FEC')
    recoveries = np.load('succinic_neutral_glucose_recoveries')
    # results_metric_4 = np.load('succinic_neutral_glucose_AOC')
    # results_metric_5 = np.load('succinic_neutral_glucose_TCI')
    inflection_product_yields = np.load('succinic_neutral_glucose_inflection_product_yields')
    
#%% succinic_neutral_sugarcane
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_sugarcane import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_neutral_sugarcane_MPSP', results_metric_1)
    np.save('succinic_neutral_sugarcane_GWP', results_metric_2)
    np.save('succinic_neutral_sugarcane_FEC', results_metric_3)
    # np.save('succinic_neutral_sugarcane_AOC', results_metric_4)
    # np.save('succinic_neutral_sugarcane_TCI', results_metric_5)
    np.save('succinic_neutral_sugarcane_recoveries', results_metric_4)
    np.save('succinic_neutral_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('succinic_neutral_sugarcane_yields', yields)
    np.save('succinic_neutral_sugarcane_titers', titers)
    np.save('succinic_neutral_sugarcane_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_neutral_sugarcane_MPSP')
    results_metric_2 = np.load('succinic_neutral_sugarcane_GWP')
    results_metric_3 = np.load('succinic_neutral_sugarcane_FEC')
    recoveries = np.load('succinic_neutral_sugarcane_recoveries')
    # results_metric_4 = np.load('succinic_neutral_sugarcane_AOC')
    # results_metric_5 = np.load('succinic_neutral_sugarcane_TCI')
    inflection_product_yields = np.load('succinic_neutral_sugarcane_inflection_product_yields')
    
#%% succinic_neutral_corn
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_corn import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_neutral_corn_MPSP', results_metric_1)
    np.save('succinic_neutral_corn_GWP', results_metric_2)
    np.save('succinic_neutral_corn_FEC', results_metric_3)
    # np.save('succinic_neutral_corn_AOC', results_metric_4)
    # np.save('succinic_neutral_corn_TCI', results_metric_5)
    np.save('succinic_neutral_corn_recoveries', results_metric_4)
    np.save('succinic_neutral_corn_inflection_product_yields', inflection_product_yields)
    np.save('succinic_neutral_corn_yields', yields)
    np.save('succinic_neutral_corn_titers', titers)
    np.save('succinic_neutral_corn_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_neutral_corn_MPSP')
    results_metric_2 = np.load('succinic_neutral_corn_GWP')
    results_metric_3 = np.load('succinic_neutral_corn_FEC')
    recoveries = np.load('succinic_neutral_corn_recoveries')
    # results_metric_4 = np.load('succinic_neutral_corn_AOC')
    # results_metric_5 = np.load('succinic_neutral_corn_TCI')
    inflection_product_yields = np.load('succinic_neutral_corn_inflection_product_yields')
    
#%% succinic_neutral_cornstover
run_simulations = True
import numpy as np
import os
product_ID = 'succinic'
if run_simulations:
    from biorefineries.succinic.analyses.TRY_analysis_cornstover import results_metric_1, results_metric_2, results_metric_3, results_metric_4,\
                                                                yields, titers, productivities,\
                                                                colors, CABBI_green_colormap,\
                                                                R302, spec, get_product_MPSP
                                                                # SA_price_range as market_range,\
    os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')
    market_range = []
    inflection_product_yields = np.array([1-0.186])
    np.save('succinic_neutral_cornstover_MPSP', results_metric_1)
    np.save('succinic_neutral_cornstover_GWP', results_metric_2)
    np.save('succinic_neutral_cornstover_FEC', results_metric_3)
    # np.save('succinic_neutral_cornstover_AOC', results_metric_4)
    # np.save('succinic_neutral_cornstover_TCI', results_metric_5)
    np.save('succinic_neutral_cornstover_recoveries', results_metric_4)
    np.save('succinic_neutral_cornstover_inflection_product_yields', inflection_product_yields)
    np.save('succinic_neutral_cornstover_yields', yields)
    np.save('succinic_neutral_cornstover_titers', titers)
    np.save('succinic_neutral_cornstover_productivities', productivities)
else:
    results_metric_1 = np.load('succinic_neutral_cornstover_MPSP')
    results_metric_2 = np.load('succinic_neutral_cornstover_GWP')
    results_metric_3 = np.load('succinic_neutral_cornstover_FEC')
    recoveries = np.load('succinic_neutral_cornstover_recoveries')
    # results_metric_4 = np.load('succinic_neutral_cornstover_AOC')
    # results_metric_5 = np.load('succinic_neutral_cornstover_TCI')
    inflection_product_yields = np.load('succinic_neutral_cornstover_inflection_product_yields')
    
    
#%% 2,3-BDO TRY
# product_ID = 'BDO'
# ## yield-titer indices should be flipped for BDO, as below
# product_ID = 'BDO'
# ind_ind_for_MPSP_y = 0
# ind_ind_for_MPSP_t = 1
# ##

# steps = 40
# yields = np.linspace(0.05, 0.95, steps)
# titers = np.linspace(5., 210., steps)
# productivities = np.array([1.00])

# market_range = []

# inflection_product_yields = [
#     1-0.055-0.02,
#     # 1-0.055,
#     ]

# os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results//BDO')

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
