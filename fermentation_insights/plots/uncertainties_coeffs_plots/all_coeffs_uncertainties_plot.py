# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:48:23 2024

@author: sarangbhagwat
"""
import numpy as np
import os
import pickle
from matplotlib import pyplot as plt
import contourplots
import itertools
from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap

#%%

#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

product_IDs = [
               'TAL', 'TAL_SA', 
               'HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol', 
               # 'succinic', 'succinic_neutral',
               ]
feedstock_IDs = [
                 'glucose', 
                 'sugarcane', 
                 'cornstover',
                 'corn', 
                 ]

all_filenames = []
coeffs_uncertainty = dict()
coeffs_baseline = dict()
refinery = {i: {} for i in feedstock_IDs}
for p,f in list(itertools.product(product_IDs, feedstock_IDs)):
    filename = p+'_'+f
    all_filenames.append(filename)
    df_dict = None
    with open(filename+'_dfd_coeffs_uncertainty.pkl', 'rb') as pickled_file:
        df_dict = pickle.load(pickled_file)
    coeffs_uncertainty[filename] = df_dict
    coeffs_baseline[filename] = cb = dict()
    cb['a'], cb['b'], cb['c'], cb['d'], cb['g'], cb['Rsq'] =\
        np.load(f'{filename}_coefficients.npy')

#%%
#%% Units
a_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}^{-1}$"
b_units = r"$\mathrm{\$}\cdot\mathrm{kg-sugars}^{-1}$"
c_units = r"$\mathrm{\$}\cdot\mathrm{L-broth}^{-1}$"
d_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}\cdot\mathrm{L-broth}^{-1}\cdot\mathrm{kg-sugars}^{-1}$"

#%%
# baseline_marker_shapes=["s", "^", "D","p", "h", "8"]
baseline_marker_shapes=["d"]*len(all_filenames)
# baseline_marker_sizes=[6, 8, 6, 10]*2
baseline_marker_sizes=[2]*len(all_filenames)
# baseline_marker_colors = ['w', '#F8858A']*4
baseline_marker_colors = ['w']*len(all_filenames)

background_fill_colors = ["#60c1cf"]*4+["#79bf82"]*4+["#90918e"]*4+["#ED586F"]*4+["#a280b9"]*4+["#f98f60"]*4
background_fill_alphas = [0.5]*24

ylabel_fontsize = 19
yticks_fontsize = 19
                         
# x_tick_labels = all_filenames 
x_tick_labels = [""]*len(all_filenames)

fig_height = 4
fig_width = 9

#%% a

a_uncertainty = [coeffs_uncertainty[filename]['a'] for filename in all_filenames]
a_baseline = [coeffs_baseline[filename]['a'] for filename in all_filenames]

contourplots.box_and_whiskers_plot(uncertainty_data=a_uncertainty, 
                          baseline_values=a_baseline,
                          baseline_marker_shapes=baseline_marker_shapes,
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_marker_colors=baseline_marker_colors,
                          baseline_locations=[i+1 for i in range(len(all_filenames))],
                          boxcolor="#A97802",
                          # ranges_for_comparison=[
                          #                        market_range,
                          #                        # [biobased_price*0.995, biobased_price*1.005],
                          #                        ],
                          # ranges_for_comparison_colors=[
                          #                               '#c0c1c2', 
                          #                               # '#646464',
                          #                               # '#c0c1c2', 
                                                        
                          #                               ],
                          # values_for_comparison=[biobased_price],
                          n_minor_ticks=3,
                          show_x_ticks=True,
                          x_tick_labels=x_tick_labels,
                          x_tick_wrap_width=9,
                          y_label=r"$\bfa$",
                          y_units=a_units,
                          # y_ticks=[f"{i:.1f}" for i in np.arange(-1, 4.01, 1)],
                          y_ticks=np.arange(-1, 4.01, 1),
                          save_file=True,
                          fig_height=fig_height,
                          fig_width=fig_width,
                          box_width=0.45,
                          filename='a_uncertainties',
                          dpi=600,
                          rotate_xticks=90.,
                          background_fill_colors=background_fill_colors,
                          background_fill_alphas=background_fill_alphas,
                          ylabel_fontsize=ylabel_fontsize,
                          yticks_fontsize=yticks_fontsize,
                          )

#%% b

# b_unc_rearr, b_base_rearr = [], []
# done_ns = []
# for f in feedstock_IDs:
#     for n in all_filenames:
#         if f in n and not f+'stover' in n and not n in done_ns:
#             print(n)
#             b_unc_rearr.append(coeffs_uncertainty[n]['b'])
#             b_base_rearr.append(coeffs_baseline[n]['b'])
#             done_ns.append(n)
            
b_uncertainty = b_unc_rearr = [coeffs_uncertainty[filename]['b'] for filename in all_filenames]
b_baseline = b_base_rearr = [coeffs_baseline[filename]['b'] for filename in all_filenames]
done_ns = all_filenames

contourplots.box_and_whiskers_plot(uncertainty_data=b_unc_rearr, 
                          baseline_values=b_base_rearr,
                          baseline_marker_shapes=baseline_marker_shapes,
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_marker_colors=baseline_marker_colors,
                          baseline_locations=[i+1 for i in range(len(all_filenames))],
                          boxcolor="#A97802",
                          # ranges_for_comparison=[
                          #                        market_range,
                          #                        # [biobased_price*0.995, biobased_price*1.005],
                          #                        ],
                          # ranges_for_comparison_colors=[
                          #                               '#c0c1c2', 
                          #                               # '#646464',
                          #                               # '#c0c1c2', 
                                                        
                          #                               ],
                          # values_for_comparison=[biobased_price],
                          n_minor_ticks=3,
                          show_x_ticks=True,
                          
                          # x_tick_labels=done_ns,
                          x_tick_labels=x_tick_labels,
                          
                          x_tick_wrap_width=9,
                          y_label=r"$\bfb$",
                          y_units=b_units,
                          # y_ticks=[f"{i:.2f}" for i in np.arange(0, 0.51, 0.1)],
                          y_ticks=np.arange(0, 0.51, 0.1),
                          save_file=True,
                          fig_height=fig_height,
                          fig_width=fig_width,
                          box_width=0.45,
                          filename='b_uncertainties',
                          dpi=600,
                          rotate_xticks=90.,
                          background_fill_colors=background_fill_colors,
                          background_fill_alphas=background_fill_alphas,
                          ylabel_fontsize=ylabel_fontsize,
                          yticks_fontsize=yticks_fontsize,)

#%% c

a_uncertainty = [coeffs_uncertainty[filename]['c'] for filename in all_filenames]
a_baseline = [coeffs_baseline[filename]['c'] for filename in all_filenames]

contourplots.box_and_whiskers_plot(uncertainty_data=a_uncertainty, 
                          baseline_values=a_baseline,
                          baseline_marker_shapes=baseline_marker_shapes,
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_marker_colors=baseline_marker_colors,
                          baseline_locations=[i+1 for i in range(len(all_filenames))],
                          boxcolor="#A97802",
                          # ranges_for_comparison=[
                          #                        market_range,
                          #                        # [biobased_price*0.995, biobased_price*1.005],
                          #                        ],
                          # ranges_for_comparison_colors=[
                          #                               '#c0c1c2', 
                          #                               # '#646464',
                          #                               # '#c0c1c2', 
                                                        
                          #                               ],
                          # values_for_comparison=[biobased_price],
                          n_minor_ticks=3,
                          show_x_ticks=True,
                          x_tick_labels=x_tick_labels,
                          x_tick_wrap_width=9,
                          y_label=r"$\bfc$",
                          y_units=c_units,
                          y_ticks=np.arange(0, 51, 10),
                          save_file=True,
                          fig_height=fig_height,
                          fig_width=fig_width,
                          box_width=0.45,
                          filename='c_uncertainties',
                          dpi=600,
                          rotate_xticks=90.,
                          background_fill_colors=background_fill_colors,
                          background_fill_alphas=background_fill_alphas,
                          ylabel_fontsize=ylabel_fontsize,
                          yticks_fontsize=yticks_fontsize,)

#%% d

a_uncertainty = [coeffs_uncertainty[filename]['d'] for filename in all_filenames]
a_baseline = [coeffs_baseline[filename]['d'] for filename in all_filenames]

contourplots.box_and_whiskers_plot(uncertainty_data=a_uncertainty, 
                          baseline_values=a_baseline,
                          baseline_marker_shapes=baseline_marker_shapes,
                          baseline_marker_sizes=baseline_marker_sizes,
                          baseline_marker_colors=baseline_marker_colors,
                          baseline_locations=[i+1 for i in range(len(all_filenames))],
                          boxcolor="#A97802",
                          # ranges_for_comparison=[
                          #                        market_range,
                          #                        # [biobased_price*0.995, biobased_price*1.005],
                          #                        ],
                          # ranges_for_comparison_colors=[
                          #                               '#c0c1c2', 
                          #                               # '#646464',
                          #                               # '#c0c1c2', 
                                                        
                          #                               ],
                          # values_for_comparison=[biobased_price],
                          n_minor_ticks=3,
                          show_x_ticks=True,
                          x_tick_labels=x_tick_labels,
                          x_tick_wrap_width=9,
                          y_label=r"$\bfd$",
                          y_units=d_units,
                          y_ticks=np.arange(-2, 2.01, 1),
                          save_file=True,
                          fig_height=fig_height,
                          fig_width=fig_width,
                          box_width=0.45,
                          filename='d_uncertainties',
                          dpi=600,
                          rotate_xticks=90.,
                          background_fill_colors=background_fill_colors,
                          background_fill_alphas=background_fill_alphas,
                          ylabel_fontsize=ylabel_fontsize,
                          yticks_fontsize=yticks_fontsize,)