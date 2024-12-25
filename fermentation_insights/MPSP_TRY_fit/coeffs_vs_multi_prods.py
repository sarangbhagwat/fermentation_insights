# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 11:48:40 2024

@author: sarangbhagwat
"""
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import contourplots
import itertools
from biosteam.utils import  colors
from  matplotlib.colors import LinearSegmentedColormap
from fermentation_insights.plots.analyze_all_combinations import coeff
from biorefineries.TAL.analyses.miscellaneous.supernatant_recycling_with_fixed_decarboxylation import plot_multiple_metrics
#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

#%%
def round_to(val, round_to):
    if isinstance(val, str):
        return val
    else:
        if round_to<1:
            return int(np.round(val,round_to))
        else:
            return np.round(val,round_to)
        
#%%
productivities = productivities_for_plot = np.linspace(0.2, 5., 7) * 0.12 # baseline productivity for TAL_SA_sugarcane = 0.12 g/L/h

product = 'TAL_SA'
feedstock = 'sugarcane'

avals, bvals, cvals, dvals, gvals = [], [], [], [], [] # note these are actually a', b', c', d', and r

for p in productivities/0.12:
    name = product+f'_{round_to(p, 1)}bp_'+feedstock if not p==1 else product+'_'+feedstock
    a, b, c, d, g, rsq = np.load(name+'_coefficients'+'.npy')
    # a = coeff[name][0]
    # b = coeff[name][1]
    # c = coeff[name][2]
    # d = coeff[name][3]
    # g = coeff[name][4]
    
    avals.append(a)
    bvals.append(b)
    cvals.append(c)
    dvals.append(d)
    gvals.append(g)

avals = np.array(avals)
bvals = np.array(bvals)
cvals = np.array(cvals)
dvals = np.array(dvals)
gvals = np.array(gvals)

#%%
##################
# Plot ticks, labels, other details
# Parameters analyzed across

x_label = r"$\bfProductivity$" # title of the x axis
# x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
x_units =r"$\mathrm{g-fp} \cdot \mathrm{L}^{-1} \cdot \mathrm{h}^{-1}$"

# Metrics
a_label = r"$\bfa'$"
b_label = r"$\bfb'$"
c_label = r"$\bfc'$"
d_label = r"$\bfd'$"
r_label = r"$\bfr$"

a_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}^{-1}$"
b_units = r"$\mathrm{\$}\cdot\mathrm{kg-sugars}^{-1}$"
c_units = r"$\mathrm{\$}\cdot\mathrm{L-broth}^{-1}$"
d_units = r"$\mathrm{\$}\cdot\mathrm{kg-fp}\cdot\mathrm{L-broth}^{-1}\cdot\mathrm{kg-sugars}^{-1}$"
r_units = r"$\mathrm{kg-p}\cdot\mathrm{kg-fp}^{-1}$"

y_labels = [a_label, b_label, c_label, d_label, r_label]
y_unitses = [a_units, b_units, c_units, d_units, r_units]

#%%
fps = 3
axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 11.
clabel_fontsize = 16
axis_tick_fontsize = 16
keep_frames = True

print('\nCreating and saving plots ...\n')

get_rounded_str = contourplots.utils.get_rounded_str

#################
 
#%%
plot_multiple_metrics(productivities,
                      [
                       avals,
                       bvals,
                       cvals,
                       dvals,
                       gvals,
                       ],
                      
                      ylabels_list=y_labels,
                      y_axis_units_list=y_unitses,
                      xlabel=x_label,
                      x_axis_units=x_units,
                      
                      units_on_newline=True,
                      xlim=(0.024,0.6), 
                      ylims_list=[
                                  (0,5),
                                  (0,1),
                                  (0,50),
                                  (-0.2, 1),
                                  (0, 0.8),
                                  ],
                      xticks=productivities,
                      yticks_list = [
                                 np.arange(0, 5.01, 1),
                                 np.arange(0, 1.01, 0.2),
                                 np.arange(0, 50.01, 10),
                                 np.arange(-0.2, 1.01, 0.2),
                                 np.arange(0, 0.801, 0.2),
                                 ],
                      y_plotline_colors_list=[
                                              "#a280b9",
                                              "#79bf82", 
                                              '#60c1cf',
                                              '#f98f60',
                                              '#90918e'
                                              ],
                      filename='coeffs_recovery_multiple_productivities_TAL_SA_sugarcane'
                      )

#%%
plot_multiple_metrics(productivities,
                      [
                       avals*productivities,
                       bvals*productivities,
                       cvals*productivities,
                       dvals*productivities,
                       gvals*productivities,
                       ],
                      
                      ylabels_list= [i + " " + r"$\bf*$" + " " + r"$\bfProductivity$" for i in y_labels],
                      y_axis_units_list=["" for i in y_unitses],
                      xlabel=x_label,
                      x_axis_units=x_units,
                      
                      units_on_newline=True,
                      xlim=(0.02,0.62), 
                      ylims_list=[
                                  (0,3),
                                  (0,0.6),
                                  (0,30),
                                  (-0.1, 0.4),
                                  (0, 0.5),
                                  ],
                      xticks=productivities,
                      yticks_list = [
                                 np.arange(0, 3.01, 0.5),
                                 np.arange(0, 0.601, 0.1),
                                 np.arange(0, 30.01, 5),
                                 np.arange(-0.1, 0.401, 0.1),
                                 np.arange(0, 0.501, 0.1),
                                 ],
                      y_plotline_colors_list=[
                                              "#a280b9",
                                              "#79bf82", 
                                              '#60c1cf',
                                              '#f98f60',
                                              '#90918e'
                                              ],
                      filename='coeffs_recovery_times_prod_multiple_productivities_TAL_SA_sugarcane'
                      )