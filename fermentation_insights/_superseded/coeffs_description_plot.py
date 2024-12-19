# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:48:40 2024

@author: sarangbhagwat
"""
import numpy as np
import os
from matplotlib import pyplot as plt
import contourplots
import itertools
from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap

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
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

#%%
def MPSP_f(y, t, a, b, c, d):
    return a + b/y + c/t + d/(y*t)

#%%
yields = yields_for_plot = np.linspace(0.05, 0.95, 50)
titers = titers_for_plot = np.linspace(10, 200, 50)

a = 0.5
b = 0.5
c = 20
d = 10.

# indicator_array = [[[MPSP_f(y, t, a, b, c, d) 
#                    for y in yields_for_plot] 
#                    for t in titers_for_plot]]

## %% 
# ########
# # corners of triangular infeasible region
# infeas_ll = 100*yields[0], titers[0]
# infeas_ul = 100*yields[0], titers[-1]

# infeas_ur_titer_index = -1
# top_titer_indicator_array = indicator_array[:, -1, :][0]
# has_infeas_region = True
# try:
#     infeas_ur_yield_index = np.where(np.isnan(top_titer_indicator_array))[0][-1] + 1
#     infeas_ur = 100*yields[infeas_ur_yield_index], titers[infeas_ur_titer_index]
# except:
#     has_infeas_region = False
    

# infeas_region_shape = {
#     # coords as tuple of tuples: (color, zorder),
#     (infeas_ll, infeas_ur, infeas_ul): ('white', 2), # infeasible region smoothing
#     } if has_infeas_region else {}

# ########

#%%
##################
# Plot ticks, labels, other details
# Parameters analyzed across

x_label = r"$\bfYield$" # title of the x axis
# x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
x_units =r"$\mathrm{g} \cdot \mathrm{g}^{-1}$"
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.]

y_label = r"$\bfTiter$" # title of the y axis
y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
y_ticks = [0, 40, 80, 120, 160, 200]

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

a_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}^{-1}$"
b_units = r"$\mathrm{\$}\cdot\mathrm{kg-sugars}^{-1}$"
c_units = r"$\mathrm{\$}\cdot\mathrm{L-broth}^{-1}$"
d_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}\cdot\mathrm{L-broth}^{-1}\cdot\mathrm{kg-sugars}^{-1}$"

#%%
fps = 3
axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 11.
clabel_fontsize = 16
axis_tick_fontsize = 16
keep_frames = True

print('\nCreating and saving contour plots ...\n')

get_rounded_str = contourplots.utils.get_rounded_str

#################

MPSP_w_levels = np.arange(0., 8.25, 0.25)
MPSP_cbar_ticks = np.arange(0., 8.1, 1.)
MPSP_w_ticks = [1., 1.25, 1.5, 1.75, 2., 2.5,  3., 4, 5, 6, 8,]
# MPSP_w_ticks = [3, 5, 8]
# MPSP_w_ticks = get_w_ticks(indicator_array, MPSP_w_levels, n_ticks=3)
# MPSP_w_ticks = get_w_ticks_from_percentiles(indicator_array, MPSP_w_levels, 
#                                             percentiles=(25, 50, 75),
#                                             cbar_ticks=MPSP_cbar_ticks)

#%% When a changes
a_vals = np.linspace(0., 5., 20)

z_label = r"$\bfCoefficient$" + " " + r"$\bfa$"# title of the z axis
z_units =  a_units
z_ticks = np.arange(a_vals[0], a_vals[-1]+1e-6, 0.5)

contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=
                                  [[[MPSP_f(y=y, t=t, a=i, b=b, c=c, d=d) 
                                   for y in yields_for_plot] 
                                   for t in titers_for_plot]
                                   for i in a_vals], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=yields_for_plot, # x axis values
                                # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                y_data=titers_for_plot, # y axis values
                                z_data=a_vals, # z axis values
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
                                fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 2),
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                extend_cmap='max',
                                cbar_ticks=MPSP_cbar_ticks,
                                z_marker_color='g', # default matplotlib color names
                                fps=fps, # animation frames (z values traversed) per second
                                n_loops=30, # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='a_coeffs_description_MPSP_y_t', # file name to save animated contourplot as (no extensions)
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
                                
                                # add_shapes = infeas_region_shape,
                                
                                units_on_newline = (False, False, False, False), # x,y,z,w
                                units_opening_brackets = [" [",] * 4,
                                units_closing_brackets = ["]",] * 4,
                                # label_over_color='white',
                                round_xticks_to=1,
                                round_yticks_to=0,
                                
                                keep_gifs=True,
                                include_top_bar = True,
                                
                                include_cbar = True,
                                include_axis_labels = True,
                                include_x_axis_ticklabels = True,
                                include_y_axis_ticklabels = True,
                                show_top_ticklabels = True,
                                figwidth=8,
                                # fig_ax_to_use=(fig, ax),
                                )

#%% When b changes
b_vals = np.linspace(0., 5., 200)

z_label = r"$\bfCoefficient$" + " " + r"$\bfb$"# title of the z axis
z_units =  b_units
z_ticks = np.arange(b_vals[0], b_vals[-1]+1e-6, 0.5)

contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=
                                  [[[MPSP_f(y=y, t=t, a=a, b=i, c=c, d=d) 
                                   for y in yields_for_plot] 
                                   for t in titers_for_plot]
                                   for i in b_vals], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=yields_for_plot, # x axis values
                                # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                y_data=titers_for_plot, # y axis values
                                z_data=b_vals, # z axis values
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
                                fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 2),
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                extend_cmap='max',
                                cbar_ticks=MPSP_cbar_ticks,
                                z_marker_color='g', # default matplotlib color names
                                fps=fps, # animation frames (z values traversed) per second
                                n_loops=30, # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='b_coeffs_description_MPSP_y_t', # file name to save animated contourplot as (no extensions)
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
                                
                                # add_shapes = infeas_region_shape,
                                
                                units_on_newline = (False, False, False, False), # x,y,z,w
                                units_opening_brackets = [" [",] * 4,
                                units_closing_brackets = ["]",] * 4,
                                # label_over_color='white',
                                round_xticks_to=1,
                                round_yticks_to=0,
                                
                                keep_gifs=True,
                                include_top_bar = True,
                                
                                include_cbar = True,
                                include_axis_labels = True,
                                include_x_axis_ticklabels = True,
                                include_y_axis_ticklabels = True,
                                show_top_ticklabels = True,
                                figwidth=8,
                                # fig_ax_to_use=(fig, ax),
                                )

#%% When c changes
c_vals = np.linspace(0., 60., 20)

z_label = r"$\bfCoefficient$" + " " + r"$\bfc$"# title of the z axis
z_units =  c_units
z_ticks = np.arange(c_vals[0], c_vals[-1]+1e-6, 5)

contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=
                                  [[[MPSP_f(y=y, t=t, a=a, b=b, c=i, d=d) 
                                   for y in yields_for_plot] 
                                   for t in titers_for_plot]
                                   for i in c_vals], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=yields_for_plot, # x axis values
                                # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                y_data=titers_for_plot, # y axis values
                                z_data=c_vals, # z axis values
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
                                fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 2),
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                extend_cmap='max',
                                cbar_ticks=MPSP_cbar_ticks,
                                z_marker_color='g', # default matplotlib color names
                                fps=fps, # animation frames (z values traversed) per second
                                n_loops=30, # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='c_coeffs_description_MPSP_y_t', # file name to save animated contourplot as (no extensions)
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
                                
                                # add_shapes = infeas_region_shape,
                                
                                units_on_newline = (False, False, False, False), # x,y,z,w
                                units_opening_brackets = [" [",] * 4,
                                units_closing_brackets = ["]",] * 4,
                                # label_over_color='white',
                                round_xticks_to=1,
                                round_yticks_to=0,
                                
                                keep_gifs=True,
                                include_top_bar = True,
                                
                                include_cbar = True,
                                include_axis_labels = True,
                                include_x_axis_ticklabels = True,
                                include_y_axis_ticklabels = True,
                                show_top_ticklabels = True,
                                figwidth=8,
                                # fig_ax_to_use=(fig, ax),
                                )

#%% When d changes
d_vals = np.linspace(0, 30., 20)

z_label = r"$\bfCoefficient$" + " " + r"$\bfd$"# title of the z axis
z_units =  d_units
z_ticks = np.arange(d_vals[0], d_vals[-1]+1e-6, 2)

contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=
                                  [[[MPSP_f(y=y, t=t, a=a, b=b, c=c, d=i) 
                                   for y in yields_for_plot] 
                                   for t in titers_for_plot]
                                   for i in d_vals], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=yields_for_plot, # x axis values
                                # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                y_data=titers_for_plot, # y axis values
                                z_data=d_vals, # z axis values
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
                                fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 2),
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                extend_cmap='max',
                                cbar_ticks=MPSP_cbar_ticks,
                                z_marker_color='g', # default matplotlib color names
                                fps=fps, # animation frames (z values traversed) per second
                                n_loops=30, # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='d_coeffs_description_MPSP_y_t', # file name to save animated contourplot as (no extensions)
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
                                
                                # add_shapes = infeas_region_shape,
                                
                                units_on_newline = (False, False, False, False), # x,y,z,w
                                units_opening_brackets = [" [",] * 4,
                                units_closing_brackets = ["]",] * 4,
                                # label_over_color='white',
                                round_xticks_to=1,
                                round_yticks_to=0,
                                
                                keep_gifs=True,
                                include_top_bar = True,
                                
                                include_cbar = True,
                                include_axis_labels = True,
                                include_x_axis_ticklabels = True,
                                include_y_axis_ticklabels = True,
                                show_top_ticklabels = True,
                                figwidth=8,
                                # fig_ax_to_use=(fig, ax),
                                )

