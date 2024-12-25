# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 17:09:07 2024

@author: sarangbhagwat
"""

#%% TAL_SA_sugarcane

#%% TAL_SA_0.2bp_sugarcane
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
#     np.save('TAL_SA_0.2bp_sugarcane_MPSP', results_metric_1)
#     np.save('TAL_SA_0.2bp_sugarcane_GWP', results_metric_2)
#     np.save('TAL_SA_0.2bp_sugarcane_FEC', results_metric_3)
#     np.save('TAL_SA_0.2bp_sugarcane_AOC', results_metric_5)
#     np.save('TAL_SA_0.2bp_sugarcane_TCI', results_metric_6)
#     np.save('TAL_SA_0.2bp_sugarcane_recoveries', results_metric_4)
#     np.save('TAL_SA_0.2bp_sugarcane_inflection_product_yields', inflection_product_yields)
#     np.save('TAL_SA_0.2bp_sugarcane_yields', yields)
#     np.save('TAL_SA_0.2bp_sugarcane_titers', titers)
#     np.save('TAL_SA_0.2bp_sugarcane_productivities', productivities)
    
# else:
#     results_metric_1 = np.load('TAL_SA_0.2bp_sugarcane_MPSP')
#     results_metric_2 = np.load('TAL_SA_0.2bp_sugarcane_GWP')
#     results_metric_3 = np.load('TAL_SA_0.2bp_sugarcane_FEC')
#     results_metric_5 = np.load('TAL_SA_0.2bp_sugarcane_AOC')
#     results_metric_6 = np.load('TAL_SA_0.2bp_sugarcane_TCI')
#     inflection_product_yields = np.load('TAL_SA_0.2bp_sugarcane_inflection_product_yields')

#%% TAL_SA_1.0bp_sugarcane
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
#     np.save('TAL_SA_1.0bp_sugarcane_MPSP', results_metric_1)
#     np.save('TAL_SA_1.0bp_sugarcane_GWP', results_metric_2)
#     np.save('TAL_SA_1.0bp_sugarcane_FEC', results_metric_3)
#     np.save('TAL_SA_1.0bp_sugarcane_AOC', results_metric_5)
#     np.save('TAL_SA_1.0bp_sugarcane_TCI', results_metric_6)
#     np.save('TAL_SA_1.0bp_sugarcane_recoveries', results_metric_4)
#     np.save('TAL_SA_1.0bp_sugarcane_inflection_product_yields', inflection_product_yields)
#     np.save('TAL_SA_1.0bp_sugarcane_yields', yields)
#     np.save('TAL_SA_1.0bp_sugarcane_titers', titers)
#     np.save('TAL_SA_1.0bp_sugarcane_productivities', productivities)
    
# else:
#     results_metric_1 = np.load('TAL_SA_1.0bp_sugarcane_MPSP')
#     results_metric_2 = np.load('TAL_SA_1.0bp_sugarcane_GWP')
#     results_metric_3 = np.load('TAL_SA_1.0bp_sugarcane_FEC')
#     results_metric_5 = np.load('TAL_SA_1.0bp_sugarcane_AOC')
#     results_metric_6 = np.load('TAL_SA_1.0bp_sugarcane_TCI')
#     inflection_product_yields = np.load('TAL_SA_1.0bp_sugarcane_inflection_product_yields')

#%% TAL_SA_1.8bp_sugarcane
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
    np.save('TAL_SA_1.8bp_sugarcane_MPSP', results_metric_1)
    np.save('TAL_SA_1.8bp_sugarcane_GWP', results_metric_2)
    np.save('TAL_SA_1.8bp_sugarcane_FEC', results_metric_3)
    np.save('TAL_SA_1.8bp_sugarcane_AOC', results_metric_5)
    np.save('TAL_SA_1.8bp_sugarcane_TCI', results_metric_6)
    np.save('TAL_SA_1.8bp_sugarcane_recoveries', results_metric_4)
    np.save('TAL_SA_1.8bp_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_1.8bp_sugarcane_yields', yields)
    np.save('TAL_SA_1.8bp_sugarcane_titers', titers)
    np.save('TAL_SA_1.8bp_sugarcane_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_1.8bp_sugarcane_MPSP')
    results_metric_2 = np.load('TAL_SA_1.8bp_sugarcane_GWP')
    results_metric_3 = np.load('TAL_SA_1.8bp_sugarcane_FEC')
    results_metric_5 = np.load('TAL_SA_1.8bp_sugarcane_AOC')
    results_metric_6 = np.load('TAL_SA_1.8bp_sugarcane_TCI')
    inflection_product_yields = np.load('TAL_SA_1.8bp_sugarcane_inflection_product_yields')

#%% TAL_SA_2.6bp_sugarcane
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
    np.save('TAL_SA_2.6bp_sugarcane_MPSP', results_metric_1)
    np.save('TAL_SA_2.6bp_sugarcane_GWP', results_metric_2)
    np.save('TAL_SA_2.6bp_sugarcane_FEC', results_metric_3)
    np.save('TAL_SA_2.6bp_sugarcane_AOC', results_metric_5)
    np.save('TAL_SA_2.6bp_sugarcane_TCI', results_metric_6)
    np.save('TAL_SA_2.6bp_sugarcane_recoveries', results_metric_4)
    np.save('TAL_SA_2.6bp_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_2.6bp_sugarcane_yields', yields)
    np.save('TAL_SA_2.6bp_sugarcane_titers', titers)
    np.save('TAL_SA_2.6bp_sugarcane_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_2.6bp_sugarcane_MPSP')
    results_metric_2 = np.load('TAL_SA_2.6bp_sugarcane_GWP')
    results_metric_3 = np.load('TAL_SA_2.6bp_sugarcane_FEC')
    results_metric_5 = np.load('TAL_SA_2.6bp_sugarcane_AOC')
    results_metric_6 = np.load('TAL_SA_2.6bp_sugarcane_TCI')
    inflection_product_yields = np.load('TAL_SA_2.6bp_sugarcane_inflection_product_yields')

#%% TAL_SA_3.4bp_sugarcane
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
    np.save('TAL_SA_3.4bp_sugarcane_MPSP', results_metric_1)
    np.save('TAL_SA_3.4bp_sugarcane_GWP', results_metric_2)
    np.save('TAL_SA_3.4bp_sugarcane_FEC', results_metric_3)
    np.save('TAL_SA_3.4bp_sugarcane_AOC', results_metric_5)
    np.save('TAL_SA_3.4bp_sugarcane_TCI', results_metric_6)
    np.save('TAL_SA_3.4bp_sugarcane_recoveries', results_metric_4)
    np.save('TAL_SA_3.4bp_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_3.4bp_sugarcane_yields', yields)
    np.save('TAL_SA_3.4bp_sugarcane_titers', titers)
    np.save('TAL_SA_3.4bp_sugarcane_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_3.4bp_sugarcane_MPSP')
    results_metric_2 = np.load('TAL_SA_3.4bp_sugarcane_GWP')
    results_metric_3 = np.load('TAL_SA_3.4bp_sugarcane_FEC')
    results_metric_5 = np.load('TAL_SA_3.4bp_sugarcane_AOC')
    results_metric_6 = np.load('TAL_SA_3.4bp_sugarcane_TCI')
    inflection_product_yields = np.load('TAL_SA_3.4bp_sugarcane_inflection_product_yields')

#%% TAL_SA_4.2bp_sugarcane
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
    np.save('TAL_SA_4.2bp_sugarcane_MPSP', results_metric_1)
    np.save('TAL_SA_4.2bp_sugarcane_GWP', results_metric_2)
    np.save('TAL_SA_4.2bp_sugarcane_FEC', results_metric_3)
    np.save('TAL_SA_4.2bp_sugarcane_AOC', results_metric_5)
    np.save('TAL_SA_4.2bp_sugarcane_TCI', results_metric_6)
    np.save('TAL_SA_4.2bp_sugarcane_recoveries', results_metric_4)
    np.save('TAL_SA_4.2bp_sugarcane_inflection_product_yields', inflection_product_yields)
    np.save('TAL_SA_4.2bp_sugarcane_yields', yields)
    np.save('TAL_SA_4.2bp_sugarcane_titers', titers)
    np.save('TAL_SA_4.2bp_sugarcane_productivities', productivities)
    
else:
    results_metric_1 = np.load('TAL_SA_4.2bp_sugarcane_MPSP')
    results_metric_2 = np.load('TAL_SA_4.2bp_sugarcane_GWP')
    results_metric_3 = np.load('TAL_SA_4.2bp_sugarcane_FEC')
    results_metric_5 = np.load('TAL_SA_4.2bp_sugarcane_AOC')
    results_metric_6 = np.load('TAL_SA_4.2bp_sugarcane_TCI')
    inflection_product_yields = np.load('TAL_SA_4.2bp_sugarcane_inflection_product_yields')

#%% TAL_SA_5.0bp_sugarcane
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
#     np.save('TAL_SA_5.0bp_sugarcane_MPSP', results_metric_1)
#     np.save('TAL_SA_5.0bp_sugarcane_GWP', results_metric_2)
#     np.save('TAL_SA_5.0bp_sugarcane_FEC', results_metric_3)
#     np.save('TAL_SA_5.0bp_sugarcane_AOC', results_metric_5)
#     np.save('TAL_SA_5.0bp_sugarcane_TCI', results_metric_6)
#     np.save('TAL_SA_5.0bp_sugarcane_recoveries', results_metric_4)
#     np.save('TAL_SA_5.0bp_sugarcane_inflection_product_yields', inflection_product_yields)
#     np.save('TAL_SA_5.0bp_sugarcane_yields', yields)
#     np.save('TAL_SA_5.0bp_sugarcane_titers', titers)
#     np.save('TAL_SA_5.0bp_sugarcane_productivities', productivities)
    
# else:
#     results_metric_1 = np.load('TAL_SA_5.0bp_sugarcane_MPSP')
#     results_metric_2 = np.load('TAL_SA_5.0bp_sugarcane_GWP')
#     results_metric_3 = np.load('TAL_SA_5.0bp_sugarcane_FEC')
#     results_metric_5 = np.load('TAL_SA_5.0bp_sugarcane_AOC')
#     results_metric_6 = np.load('TAL_SA_5.0bp_sugarcane_TCI')
#     inflection_product_yields = np.load('TAL_SA_5.0bp_sugarcane_inflection_product_yields')
