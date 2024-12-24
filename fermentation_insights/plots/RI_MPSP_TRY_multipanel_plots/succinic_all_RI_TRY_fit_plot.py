# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 16:48:23 2024

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

def GG_orange_white_blue_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.

    """
    # CABBI_colors = (colors.CABBI_orange.RGBn,
    #                 colors.CABBI_yellow.RGBn,

    #                 colors.CABBI_green.RGBn,
    #                 # colors.CABBI_teal_green.shade(50).RGBn,
    #                 colors.grey_dark.RGBn)
    cmap_colors = np.array((
                   (99, 198, 206), # GG blue
                   (138, 227, 235),# lighter GG blue
                   
                   (255, 255, 255), # white
                   
                   (250, 161, 82), # lighter GG orange
                   (229, 135, 53), # GG orange
                   ))/255
    
    return LinearSegmentedColormap.from_list('CABBI', cmap_colors, N_levels)

#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

product_IDs = [
               # 'TAL', 'TAL_SA', 
               # 'HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol', 
                'succinic', 'succinic_neutral',
               ]
feedstock_IDs = [
                 'glucose', 
                 'cornstover',
                 'sugarcane', 
                 'corn', 
                 ]

all_filenames = []
refinery = {i: {} for i in feedstock_IDs}
for p,f in list(itertools.product(product_IDs, feedstock_IDs)):
    all_filenames.append(p+'_'+f)

#%%
fig, axs = plt.subplots(len(product_IDs), len(feedstock_IDs), constrained_layout=True)

#%% Plot clabel utils
def get_penultimate_lowest_w_tick(w_levels, w_array):
    return_next_i = False
    for i in w_levels:
        if return_next_i: 
            return i
        elif np.any(w_array<i):
            return_next_i = True
    raise ValueError('\nw_array does not contain any values lower than the second-highest value in cbar_ticks.\n')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def get_w_ticks_from_bounds(w_array, w_levels, n_ticks=3):   # assumes w_array contains values lower than the second-highest value in w_levels
    low_bound = get_penultimate_lowest_w_tick(w_levels, w_array)
    high_bound = w_levels[-1]
    # w_ticks_unrefined = np.linspace(low_bound, high_bound+np.min(w_levels)/10, n_ticks)
    w_ticks_unrefined = np.linspace(low_bound, high_bound, n_ticks)
    w_ticks_refined = []
    for i in w_ticks_unrefined:
        new_w_tick = find_nearest(w_levels, i)
        if new_w_tick not in w_ticks_refined:
            w_ticks_refined.append(new_w_tick)
    return w_ticks_refined

def get_w_ticks_from_percentiles(w_array, w_levels, percentiles=(25, 50, 75), cbar_ticks=None):   # assumes w_array contains values lower than the second-highest value in w_levels
    w_arr_non_nan = w_array[np.where(~np.isnan(w_array))]
    if cbar_ticks is not None:
        w_arr_non_nan = w_arr_non_nan[np.where(w_arr_non_nan<np.max(cbar_ticks))]
    w_ticks_refined = []
    for p in percentiles:
        try:
            i = np.percentile(w_arr_non_nan, p)
        except:
            breakpoint()
        new_w_tick = find_nearest(w_levels, i)
        if new_w_tick not in w_ticks_refined:
            w_ticks_refined.append(new_w_tick)
    return w_ticks_refined

#%% Plot

log_scale_ri = True

for i, product_ID in zip(range(len(product_IDs)), product_IDs):
    axs_list = None
    if len(product_IDs)>1:
        axs_list= axs[i]
    else:
        axs_list = [axs[i]]
        
    for j, feedstock_ID in zip(range(len(feedstock_IDs)), feedstock_IDs):
        ax = None
        if len(feedstock_IDs)>1:
            ax = axs_list[j]
        else:
            ax = axs_list
        filename = product_ID+'_'+feedstock_ID
        print(filename)
        yields = yields_for_plot = np.load(filename+'_yields_for_eval.npy')
        titers = titers_for_plot = np.load(filename+'_titers_for_eval.npy')
        RI_array = results_metric_1 = np.load(f'{filename}_RI_MPSP.npy')
        if log_scale_ri: RI_array = np.log10(RI_array)
        
        show_yields = product_ID in ['TAL_SA', 'HP_neutral_hexanol', 'succinic_neutral']
        show_titers = feedstock_ID in['glucose']
        show_top_ticklabels = ((product_ID in ['TAL', 'HP', 'succinic'] and feedstock_ID in ['glucose'])
                               or product_ID in ['TAL_SA', 'HP_neutral_hexanol', 'succinic_neutral'] and feedstock_ID in ['cornstover'])
        
        ########
        # corners of triangular infeasible region
        infeas_ll = 100*yields[0], titers[0]
        infeas_ul = 100*yields[0], titers[-1]
        
        infeas_ur_titer_index = -1
        top_titer_RI_array = RI_array[:, -1, :][0]
        has_infeas_region = True
        try:
            infeas_ur_yield_index = np.where(np.isnan(top_titer_RI_array))[0][-1] + 1
            infeas_ur = 100*yields[infeas_ur_yield_index], titers[infeas_ur_titer_index]
        except:
            has_infeas_region = False
            
        
        infeas_region_shape = {
            # coords as tuple of tuples: (color, zorder),
            (infeas_ll, infeas_ur, infeas_ul): ('white', 2), # infeasible region smoothing
            } if has_infeas_region else {}
        
        ########
        
        ##################
        # Plot ticks, labels, other details
        # Parameters analyzed across

        x_label = r"$\bfYield$" # title of the x axis
        # x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
        x_units =r"$\mathrm{g} \cdot \mathrm{g}^{-1}$"
        # x_ticks = [20, 30, 40, 50, 60, 70, 80, 90, 100]

        y_label = r"$\bfTiter$" # title of the y axis
        y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"

        # ##### X and Y ticks -- bounds rounding method (use for yield)        
        round_up_yield_ubound = True if np.round(yields_for_plot[-1],1)<yields_for_plot[-1] else False
        round_down_yield_lbound = True if np.round(yields_for_plot[0],1)>yields_for_plot[0] else False
        x_ticks = np.arange(np.round(yields_for_plot[0],1)-0.1*round_down_yield_lbound, 
                            np.round(yields_for_plot[-1],1)+0.1*round_up_yield_ubound+1e-5, 0.1)[::4]


        # round_up_titer_ubound = True if np.round(titers_for_plot[-1],-1)<titers_for_plot[-1] else False
        # # round_down_titer_lbound = True if np.round(titers_for_plot[0],-1)>titers_for_plot[0] else False
        # y_ticks = np.arange(np.round(0,-1), 
        #                     np.round(titers_for_plot[-1],-1)+10*round_up_titer_ubound+1e-5, 10)[::4]
        # #####
        
        
        ##### X and Y ticks -- bounds and remainder method (use for titer)
        titer_ub = titers_for_plot[-1]
        titer_lb = 0.
        n_steps = None
        if titer_ub%4==0:
            n_steps = 5
        elif titer_ub%3==0:
            n_steps = 4
        elif titer_ub%2==0:
            n_steps = 3
        y_ticks = np.linspace(titer_lb, titer_ub, n_steps)
        
        # yield_ub = yields_for_plot[-1]
        # yield_lb = 0.
        # n_steps = None
        # if yield_ub*100%4==0:
        #     n_steps = 5
        # elif yield_ub*100%3==0:
        #     n_steps = 4
        # elif yield_ub*100%2==0:
        #     n_steps = 3
        # x_ticks = np.linspace(yield_lb, yield_ub, n_steps)
        
        #####
        
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
        clabel_fontsize = 16
        axis_tick_fontsize = 16
        keep_frames = True

        print('\nCreating and saving contour plots ...\n')

        get_rounded_str = contourplots.utils.get_rounded_str
        
        #################
        
        rel_impact_w_levels = np.arange(-3, 3.01, 0.25)
        rel_impact_cbar_ticks = np.arange(-3, 3.01, 0.5)
        # rel_impact_w_ticks = [1., 1.25, 1.5, 1.75, 2., 2.5,  3., 4, 5, 6, 8,]
        rel_impact_w_ticks = [-1, 0, 1, 1.5, 2]
        # rel_impact_w_ticks = get_w_ticks(RI_array, rel_impact_w_levels, n_ticks=3)
        # rel_impact_w_ticks = get_w_ticks_from_percentiles(RI_array, rel_impact_w_levels, 
        #                                             percentiles=(25, 50, 75),
        #                                             # cbar_ticks=rel_impact_cbar_ticks,
        #                                             cbar_ticks=None,
        #                                             )
        
        contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=RI_array, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., rel_impact
                                        x_data=yields_for_plot, # x axis values
                                        # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                        y_data=titers_for_plot, # y axis values
                                        z_data=np.array([1]), # z axis values
                                        x_label=x_label, # title of the x axis
                                        y_label=y_label, # title of the y axis
                                        z_label=z_label, # title of the z axis
                                        w_label=rel_impact_w_label, # title of the color axis
                                        x_ticks=x_ticks,
                                        y_ticks=y_ticks,
                                        z_ticks=z_ticks,
                                        w_levels=rel_impact_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                        w_ticks=rel_impact_w_ticks, # labeled, lined contours; a subset of w_levels
                                        x_units=x_units,
                                        y_units=y_units,
                                        z_units=z_units,
                                        w_units=rel_impact_units,
                                        # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                        fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 2),
                                        cmap=GG_orange_white_blue_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                        cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                        extend_cmap='neither',
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
                                        
                                        add_shapes = infeas_region_shape,
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" (",] * 4,
                                        units_closing_brackets = [")",] * 4,
                                        # label_over_color='white',
                                        round_xticks_to=1,
                                        round_yticks_to=0,
                                        keep_gifs=False,
                                        include_top_bar = False,
                                        include_cbar = False,
                                        include_axis_labels = False,
                                        include_x_axis_ticklabels = show_yields,
                                        include_y_axis_ticklabels = show_titers,
                                        show_top_ticklabels = show_top_ticklabels,
                                        figwidth=None,
                                        fig_ax_to_use=(fig, ax),
                                        inline_spacing=0.,
                                        )
        
        # # add R-squared text box to top left of each plot
        # coefficients = np.load(f'{filename}_coefficients.npy')
        # Rsq = coefficients[-1]
        # textstr = "$\mathrm{R}^{2}$" + " = " + "%.3f"%(Rsq)
        # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        
        # # place a text box in upper left in axes coords
        # ax.text(
        #         np.linspace(x_ticks[0], x_ticks[-1], 100)[4], # x-coord
        #         # np.linspace(x_ticks[0], x_ticks[-1], 100)[15], # x-coord
        #         np.linspace(y_ticks[0], y_ticks[-1], 100)[-5], # y-coord
        #         textstr, 
        #         # transform=ax.transAxes, 
        #         fontsize=axis_tick_fontsize,
        #         # verticalalignment='top', 
        #         bbox=props)
        
#%%
plt.subplots_adjust(wspace=0, hspace=0)
fig.set_figheight(8.5)
fig.set_figheight(4.5)
fig.set_figwidth(9)

plt.savefig(f'succinic_all_RI_MPSP_TRY_fit.png', 
            transparent = False,  
            facecolor = 'white',
            bbox_inches='tight',
            dpi=300,
            )  
