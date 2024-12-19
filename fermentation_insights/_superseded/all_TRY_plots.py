# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 15:28:12 2024

@author: sarangbhagwat
"""

import numpy as np
import contourplots
from matplotlib import pyplot as plt

#%%

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
MPSP_w_ticks = [1., 1.5, 1.75, 2., 2.5,  3., 3.5, 4, 6, 8,]
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