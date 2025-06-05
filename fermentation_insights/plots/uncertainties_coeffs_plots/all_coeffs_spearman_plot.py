# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 19:29:00 2025

@author: sarangbhagwat
"""

import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.collections import EllipseCollection
from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap

os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

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
#%% Utils
def plot_corr_ellipses(data, 
                       ax=None,  
                       size_factor=0.85,
                       circles=True,
                       threshold_abs_value=0.1,
                       **kwargs,):

    M = list(np.array(data))
    # # add rows of ones and negative ones to get color levels -1 to 1; cannot set levels for collections
    # M.append(np.ones(len(M[-1])))
    # M.append(-1*np.ones(len(M[-1])))
    
    M = np.array(M)
    if not M.ndim == 2:
        raise ValueError('data must be a 2D array')
    if ax is None:
        fig, ax = plt.subplots(1, 1, subplot_kw={'aspect':'equal'})
        ax.set_xlim(-0.5, M.shape[1] - 0.5)
        ax.set_ylim(-0.5, M.shape[0] - 0.5)
    
    ax.set_xticks(np.arange(-0.5, M.shape[1] - 0.5 + 1e-6, 1.))
    ax.set_yticks(np.arange(-0.5, M.shape[0] - 0.5 + 1e-6, 1.))
    
    # xy locations of each ellipse center
    xy = np.indices(M.shape)[::-1].reshape(2, -1).T

    # set the relative sizes of the major/minor axes according to the strength of
    # the positive/negative correlation
    if circles:
        w = np.abs(M).ravel()
        h = np.abs(M).ravel()
    else:
        w = np.ones_like(M).ravel()
        h = 1 - np.abs(M).ravel()
    
    if threshold_abs_value:
        w[np.where(np.abs(w)<threshold_abs_value)] = 0.
        h[np.where(np.abs(h)<threshold_abs_value)] = 0.
        
    a = 45 * np.sign(M).ravel()

    ec = EllipseCollection(widths=w*size_factor, heights=h*size_factor, 
                           angles=a, units='x', offsets=xy,
                           transOffset=ax.transData, array=M.ravel(), **kwargs)
    ax.add_collection(ec)
    
    plt.grid(axis='both', which='major',
             # color='0.95',
             )
    
    # if data is a DataFrame, use the row/column names as tick labels
    if isinstance(data, pd.DataFrame):
        ax.set_xticks(np.arange(M.shape[1]))
        ax.set_xticklabels(data.columns, rotation=90)
        ax.set_yticks(np.arange(M.shape[0]))
        ax.set_yticklabels(data.index)
    
    return ec

#%% Get data
takeaways = {}
takeaways_by_param = {}

params_in_desired_order = [
 'operating time',
 'feedstock unit price',
 'feedstock capacity',
 'inoculum ratio',
 'seed train fermentation ratio',
 'fermentation cell mass yield',
 'CSL unit price',
 'CSL loading',
 'product storage time',
 'natural gas unit price',
 'boiler efficiency',
 'electricity unit price',
 'turbogenerator efficiency',]

for coeff in ['a', 'b', 'c', 'd']:
    takeaways_by_param[coeff] = []
    takeaways[coeff] = pd.read_excel(coeff+'_spearman_general_takeaways.xlsx')
    for param in ['Unnamed: 0', 'product', 'feedstock',]:
        takeaways[coeff] = takeaways[coeff].drop(param, axis=1)
    # for param in takeaways[coeff].columns:
    params_in_desired_order_reversed = params_in_desired_order.copy()
    params_in_desired_order_reversed.reverse()
    for param in params_in_desired_order_reversed:
            takeaways_by_param[coeff].append(takeaways[coeff][param])

#%% Plot

coeff = 'b'

plt.rcParams['font.sans-serif'] = "Arial Unicode"
plt.rcParams['font.size'] = str(20)

fig, ax = plt.subplots(1, 1, constrained_layout=True)

data = takeaways_by_param[coeff]


m = plot_corr_ellipses(data, ax=ax, 
                       # cmap='Greens',
                       cmap=GG_orange_white_blue_colormap(),
                       clim=[-1,1]
                       )

ax.margins(0.1)
cb = fig.colorbar(m, ax=ax,
                  ticks=[-1.0, -0.8, -0.6, -0.4, -0.2,
                         0,
                         0.2, 0.4, 0.6, 0.8, 1.0],
                  format='%.1f'
                  )
cb.ax.minorticks_on()

# cb = plt.colorbar(m, ax=ax, ticks = [-1, -0.5, 0, 0.5, 1.])
cb.set_label('Correlation coefficient')

# Save figure
fig.set_figheight(12)
fig.set_figwidth(24.)

plt.savefig(f'{coeff}_spearman_takeaways.png', 
            transparent=False,  
            facecolor='white',
            bbox_inches='tight',
            dpi=300,
            )  
